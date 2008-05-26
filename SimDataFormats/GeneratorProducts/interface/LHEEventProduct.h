#ifndef SimDataFormats_GeneratorProducts_LHEEventProduct_h
#define SimDataFormats_GeneratorProducts_LHEEventProduct_h

#include <memory>
#include <vector>
#include <string>

#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"

class LHEEventProduct {
    public:
	struct PDF {
		std::pair<int, int>		id;
		std::pair<double, double>	x;
		std::pair<double, double>	xPDF;
		double				scalePDF;
	};

	typedef std::vector<std::string>::const_iterator
						comments_const_iterator;
	typedef std::vector<std::string>::size_type size_type;

	LHEEventProduct() {}
	LHEEventProduct(const lhef::HEPEUP &hepeup) : hepeup_(hepeup) {}
	~LHEEventProduct() {}

	void setPDF(const PDF &pdf) { pdf_.reset(new PDF(pdf)); }
	void addComment(const std::string &line) { comments_.push_back(line); }

	const lhef::HEPEUP &hepeup() const { return hepeup_; }
	const PDF *pdf() const { return pdf_.get(); }

	size_type comments_size() const { return comments_.size(); }
	comments_const_iterator comments_begin() const { return comments_.begin(); }
	comments_const_iterator comments_end() const { return comments_.end(); }

	class const_iterator {
	    public:
		typedef std::forward_iterator_tag	iterator_category;
		typedef std::string			value_type;
		typedef std::ptrdiff_t			difference_type;
		typedef std::string			*pointer;
		typedef std::string			&reference;

		const_iterator() : line(npos) {}
		~const_iterator() {}

		inline bool operator == (const const_iterator &other) const
		{ return line == other.line; }
		inline bool operator != (const const_iterator &other) const
		{ return !operator == (other); }

		inline const_iterator &operator ++ ()
		{ next(); return *this; }
		inline const_iterator operator ++ (int dummy)
		{ const_iterator orig = *this; next(); return orig; }

		const std::string &operator * () const { return tmp; }
		const std::string *operator -> () const { return &tmp; }

	    private:
		friend class LHEEventProduct;

		void next();

		const LHEEventProduct	*event;
		unsigned int		line;
		std::string		tmp;

		static const unsigned int npos = 99999;
	};

	inline const_iterator begin() const
	{ const_iterator iter; iter.event = this; iter.line = 0; return iter; }
	inline const_iterator end() const { return const_iterator(); }

    private:
	lhef::HEPEUP			hepeup_;
	std::vector<std::string>	comments_;
	std::auto_ptr<PDF>		pdf_;
};

#endif // GeneratorEvent_LHEInterface_LHEEventProduct_h
