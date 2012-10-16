/*
 * SharedQueue.h
 */

#ifndef _SHARED_QUEUE_H
#define _SHARED_QUEUE_H

#include "Lock.h"
#include "compiler.h"
#include <stddef.h>
#include <string.h>

//#define SHARED_QUEUE_LOCKFREE

/*
 * A shared producer-consumer queue implementation with blocking put()s and
 * get()s, and size determined at construction.
 *
 * The type T must be a pointer type, and NULL's cannot be added to the queue.
 *
 * If SHARED_QUEUE_LOCKFREE is defined, the put() and get() methods will use
 * compare-and-swap operations to avoid locking @m_lock.  Threads may still have
 * to sleep to wait on the appropriate semaphore; however, there will be less
 * thread contention due to not having to acquire @m_lock as well.
 *
 * For the lock-free queue, it cannot be guaranteed that the pointers placed in
 * the m_array are actually between m_front and m_back.  For example, the queue
 * could wrap all the way around the array such that when a get()ing thread
 * wakes up, it picks up something from the middle of the valid pointers, thus
 * creating a hole in the queue.  However, the code should not break when this
 * happens, and repeated put()'s should usually place the pointers contiguously,
 * allowing the holes to eventually be gotten rid of once the queue wraps around
 * again.
 *
 * It is also possible for the ordering to be violated in the lock-free version
 * (i.e. in a very unlikely situation, consecutive get()s may return values in
 * an order not the same as the order they were inserted.).  However, this can
 * only happen if multiple threads are get()ing or put()ing values concurrently.
 */
template <typename T>
class SharedQueue {
private:
	Semaphore  m_filled_slots;
	Semaphore  m_empty_slots;
	unsigned   m_front;
	unsigned   m_back;
	T	  *m_array;
#ifndef SHARED_QUEUE_LOCKFREE
	Lock	   m_lock;
#endif
	size_t	   m_size;

#ifdef SHARED_QUEUE_LOCKFREE

	/*
	 * Atomically increment the unsigned integer pointed to by @var,
	 * wrapping to 0 when the queue size is reached.  @old_val specifies
	 * what we think the current value of the variable is (since we read it
	 * earlier, we don't need to read it again until the actual CAS
	 * operation).
	 */
	inline void increment_bounds_variable(volatile unsigned *var,
					      unsigned old_val)
	{
		size_t size = m_size;
		unsigned current_val, new_val;
		while (1) {
			new_val = (old_val + 1) % size;
			current_val = cas(var, old_val, new_val);
			if (current_val == old_val)
				return;
			old_val = current_val;
		}
	}
#endif

public:

	SharedQueue(size_t size)
		: m_filled_slots(0),
		  m_empty_slots(size),
		  m_front(0),
		  m_back(size - 1),
		  m_array(new T[size]),
#ifndef SHARED_QUEUE_LOCKFREE
		  m_lock(),
#endif
		  m_size(size)
	{ memset(m_array, 0, size * sizeof(T)); }

	~SharedQueue() {
		for (size_t i = 0; i < m_size; i++)
			delete m_array[i];
		delete [] m_array;
	}


#ifdef SHARED_QUEUE_LOCKFREE
	/* Retrieve a T from the queue as soon as there is one available. */
	T get() {
		T t;
		unsigned front, old_front;
		size_t size = m_size;

		m_filled_slots.wait();
		old_front = front = *(volatile unsigned*)&m_front;
		while (1) {
			t = m_array[front]; // Pointer we are going to try to take
			if (t != NULL && cas_bool(&m_array[front], t, (T)NULL))
				break; // We claimed the pointer.
			// Someone else claimed the pointer first, or there is
			// no valid pointer here.  Try the next slot.
			front = (front + 1) % size;
		}

		// Advance m_front to the next slot (this sometimes may not
		// actually be the slot following the pointer we took)
		increment_bounds_variable((volatile unsigned *)&m_front, old_front);

		m_empty_slots.post();
		return t;
	}
	/* Place the T in the queue as soon as there is an empty space. */
	void put(T t) {
		unsigned back, old_back;
		size_t size = m_size;

		m_empty_slots.wait();
		
		old_back = back = *(volatile unsigned*)&m_back;
		back = (back + 1) % size;
		while (!cas_bool(&m_array[back], (T)NULL, t)) {
			// Someone else put() a T into the slot first.  Try the
			// next one.
			back = (back + 1) % size;
		}
		// Advance m_back to the next slot (this sometimes may not
		// actually be the slot corresponding to the pointer we placed
		// in the array)
		increment_bounds_variable((volatile unsigned *)&m_back, old_back);

		m_filled_slots.post();
	}
#else // SHARED_QUEUE_LOCKFREE

	/* Retrieve a T from the queue as soon as there is one available. */
	T get() {
		m_filled_slots.wait();
		m_lock.lock();

		T t = m_array[m_front];
		m_array[m_front] = (T)NULL;
		m_front = (m_front + 1) % m_size;

		m_empty_slots.post();
		m_lock.unlock();
		return t;
	}
	/* Place the T in the queue as soon as there is an empty space. */
	void put(T t) {
		m_empty_slots.wait();
		m_lock.lock();

		m_back = (m_back + 1) % m_size;
		m_array[m_back] = t;

		m_filled_slots.post();
		m_lock.unlock();
	}
#endif // SHARED_QUEUE_LOCKFREE
};

#endif // _SHARED_QUEUE_H
