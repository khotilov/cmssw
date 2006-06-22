#ifndef CSCSPTrailer_h
#define CSCSPTrailer_h

class CSCSPTrailer {
private:
	/////// word 1 ///////
	unsigned l1a_           : 8;
	unsigned word_count_low : 4;
	unsigned trailer_mark_1 : 4;  // constant, should be 1111 = 0xF
	/////// word 2 ///////
	unsigned trailer_mark_2 : 4;  // constant, should be 1111 = 0xF
	unsigned trailer_mark_3 : 3;  // constant, should be  111 = 0x7
	unsigned l1a_fifo_full  : 1;
	unsigned word_count_high: 4;
	unsigned trailer_mark_4 : 4;  // constant, should be 1111 = 0xF
	/////// word 3 ///////
	unsigned core_release_month: 4;
	unsigned core_release_year : 4;
	unsigned spare_1           : 1;
	unsigned spare_2           : 1;
	unsigned spare_3           : 1;
	unsigned zero_1            : 1;
	unsigned trailer_mark_5    : 4;  // constant, should be 1111 = 0xF
	/////// word 4 ///////
	unsigned core_configuraton : 12;
	unsigned trailer_mark_6    : 4;  // constant, should be 1111 = 0xF
	/////// word 5 ///////
	unsigned core_release_day  : 5;
	unsigned zero_2            : 7;
	unsigned trailer_mark_7    : 4;  // constant, should be 1110 = 0xE
	/////// word 6 ///////
	unsigned board_id          : 12;
	unsigned trailer_mark_8    : 4;  // constant, should be 1110 = 0xE
	/////// word 7 ///////
	unsigned crc_low           : 11;
	unsigned crc_low_parity    : 1;
	unsigned trailer_mark_9    : 4;  // constant, should be 1110 = 0xE
	/////// word 8 ///////
	unsigned crc_high          : 11;
	unsigned crc_high_parity   : 1;
	unsigned trailer_mark_10   : 4;  // constant, should be 1110 = 0xE

public:
	bool check(void) const throw() {
		return spare_1!=0 || spare_2!=0 || spare_3!=0 || zero_1!=0 || zero_2!=0 ||
			trailer_mark_1!=0xF || trailer_mark_2!=0xF || trailer_mark_3!=0x7 || trailer_mark_4!=0xF || trailer_mark_5!=0xF || trailer_mark_6!=0xF ||
			trailer_mark_7!=0xE || trailer_mark_8!=0xE || trailer_mark_9!=0xE || trailer_mark_10!=0xE;
	}

	bool unpack(const unsigned short *&buf) throw()  { memcpy(this, buf, 8*sizeof(short)); buf+=8; return check(); }

	CSCSPTrailer(void){}
};

#endif
