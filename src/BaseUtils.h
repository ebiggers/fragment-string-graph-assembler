#pragma once

class BaseUtils {
private:
	static const unsigned char _ascii_to_bin_tab[256];
	static const unsigned char _bin_to_ascii_tab[4];
public:

	static inline unsigned char ascii_to_bin(char ascii_base) {
		return _ascii_to_bin_tab[static_cast<unsigned char>(ascii_base)];
	}
	static inline char bin_to_ascii(unsigned char bin_base) {
		return _bin_to_ascii_tab[bin_base];
	}
};
