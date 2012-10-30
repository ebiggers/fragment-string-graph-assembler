#pragma once

class BaseUtils {
private:
	static const unsigned char _ascii_to_bin_tab[256];
public:

	static inline unsigned char ascii_to_bin(char ascii_base) {
		return _ascii_to_bin_tab[static_cast<unsigned char>(ascii_base)];
	}
};
