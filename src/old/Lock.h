/*
 * Lock.h
 */

#ifndef _LOCK_H
#define _LOCK_H

#include <pthread.h>
#include <semaphore.h>

class Lock {
private:
	pthread_mutex_t _lock;
public:
	Lock()	     { pthread_mutex_init(&_lock, NULL);}
	~Lock()	     { pthread_mutex_destroy(&_lock);	}
	void set()   { pthread_mutex_lock(&_lock);	}
	void lock()  { set();				}
	void unset() { pthread_mutex_unlock(&_lock);	}
	void unlock(){ unset();				}
	bool test()  { return pthread_mutex_trylock(&_lock) == 0; }
};

class Cond
{
private:
	pthread_mutex_t  mutex;
	pthread_cond_t   cond;

public:
	Cond() { 
		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);
	}

	~Cond() {
		pthread_cond_destroy(&cond);
		pthread_mutex_destroy(&mutex);
	}

	void lock()	 { pthread_mutex_lock(&mutex);		}
	void unlock()	 { pthread_mutex_unlock(&mutex);	}
	void wait()	 { pthread_cond_wait(&cond, &mutex);	}
	void signal()	 { pthread_cond_signal(&cond);		}
	void broadcast() { pthread_cond_broadcast(&cond);	}
};

class Semaphore {
private:
	sem_t sem;
public:
	Semaphore() {
		sem_init(&sem, 0, 1);
	}

	Semaphore(unsigned val) {
		sem_init(&sem, 0, val);
	}

	~Semaphore() {
		sem_destroy(&sem);
	}

	void post() {
		sem_post(&sem);
	}

	void wait() {
		sem_wait(&sem);
	}

	int val() {
		int n;
		sem_getvalue(&sem, &n);
		return n;
	}
	operator int() {
		return val();
	}
};

#endif
