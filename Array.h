#ifndef ARRAY_H
#define ARRAY_H
#include <iostream>
#include <initializer_list>
#include <cassert>
#include <complex>
#include <iomanip>
#include <cstddef>
#include "mkl.h"
#include <iterator>
#include <string>
#include <fstream>
#include <limits>
#include <type_traits>
#include <random>
#include <utility>
#include <functional>
// namespace name: act (array class template)
// A array class template

// Add a template specialization for to_file functions when T == std::complex<double, float>
namespace std {
	template <>
	class numeric_limits<std::complex<double>> : public numeric_limits<double> {};
	template <>
	class numeric_limits<std::complex<float>> : public numeric_limits<float> {};
}

namespace act {
	// Alias for complex type.
	template<typename T>
	using cx = std::complex<T>;

	//inline unsigned int random_seed{rd()};
	//unsigned int random_seed;

	template<typename T>
	struct Uniform_distribution {
		static auto uniform_distribution(T from, T to);
	};

	template<typename T>
	auto Uniform_distribution<T>::uniform_distribution(T from, T to) {
		return std::uniform_real_distribution<T>(from, to);
	}

	template<>
	struct Uniform_distribution<int> {
		static auto uniform_distribution(int from, int to);
	};
	auto Uniform_distribution<int>::uniform_distribution(int from, int to) {
		return std::uniform_int_distribution<>(from, to);
	}

	// Return a function which takes a mersenne and gives a complex.
	template<typename T>
	std::function<act::cx<T>(std::mt19937_64&)> uniform_complex_distribution(act::cx<T> from, act::cx<T> to) {
		std::uniform_real_distribution<T> dr(std::real(from), std::real(to));
		std::uniform_real_distribution<T> di(std::imag(from), std::imag(to));

		auto res = 
			[dr, di]
			(std::mt19937_64& n_gen) mutable 
		{
			act::cx<T> value;
			value.real(dr(n_gen)); 
			value.imag(di(n_gen));
			return value;
		};

		return res;
	}

	template<typename T>
	struct Uniform_distribution<act::cx<T>> {
		static auto uniform_distribution(act::cx<T> from, act::cx<T> to);
	};

	template<typename T>
	auto act::Uniform_distribution<act::cx<T>>::uniform_distribution(act::cx<T> from, act::cx<T> to) {
		return uniform_complex_distribution(from, to);
	}
	

	template<typename T>
	class Matrix;

	template<typename T>
	class Array;


	// This function initializes the random generator for a given type.
	// If not used, the generator is randomly initialized.
	template<typename T>
	static void set_seed(unsigned int seed);

	template<typename S>
	Array<S> operator*(Matrix<S> const& M, Array<S> const& A);

	template<typename T>
	class Array {
	protected:
		T * ptr;
		static std::mt19937_64 r_generator;
		std::size_t size;
	public:
		friend class act::Matrix<T>;

		// Iterators implementation.

		class iterator {
		private:
			friend class Array<T>;
			using value_type = T;
			using difference_type = std::ptrdiff_t;
			using pointer = T *;
			using reference = T&;
			using iterator_category = std::random_access_iterator_tag;
			pointer ptr_;
			std::size_t index_;
		public:
			iterator() : ptr_(nullptr), index_(0) {}
			iterator(pointer ar, std::size_t index) : ptr_(ar), index_(index) {}
			iterator(iterator const& itr) : ptr_(itr.ptr_), index_(itr.index_) {}

			~iterator() {}

			iterator& operator=(iterator const& itr) {
				if (this == &itr)
					return *this;
				this->ptr_ = itr.ptr_;
				this->index_ = itr.index_;
				return *this;
			}

			reference operator*() {
				return this->ptr_[index_];
			}
			
			pointer operator->() {
				return this->ptr_;
			}

			iterator& operator++() {
				++this->index_;
				return *this;
			}
			iterator operator++(int) {
				iterator itr(this->ar, this->index_);
				this->index_++;
				return itr;
			}
			iterator& operator--() {
				this->index_--;
				return *this;
			}
			iterator operator--(int) {
				iterator itr(this->ar, this->index_);
				this->index_--;
				return itr;
			}
			friend iterator operator+(iterator const& itr1, iterator const& itr2) {
				assert(&itr1.ar == &itr2.ar);
				iterator itr3(itr1.ar, itr1.index_ + itr2.index_);
				return itr3;
			}
			friend iterator operator-(iterator const& itr1, iterator const& itr2) {
				assert(&itr1.ar == &itr2.ar);
				iterator itr3(itr1.ar, itr1.index_ - itr2.index_);
				return itr3;
			}
			friend iterator operator+(iterator const& itr1, difference_type num) {
				iterator itr2(itr1.ar, itr1.index_ + num);
				return itr2;
			}
			friend iterator operator+(difference_type num, iterator const& itr1) {
				return itr1 + num;
			}
			friend iterator operator-(difference_type num, iterator const& itr1) {
				iterator itr2(itr1.ptr_, itr1.index_ - num);
				return itr2;
			}
			friend iterator operator-(iterator const& itr1, difference_type num) {
				iterator itr2(itr1.ptr_, itr1.index_ - num);
				return itr2;
			}
			iterator& operator+=(difference_type num) {
				this->index_ += num;
				return *this;
			}
			iterator& operator-=(difference_type num) {
				this->index_ -= num;
				return *this;
			}
			reference operator[](difference_type ind) {
				return this->ptr_[ind];
			}
			friend void swap(iterator itr1, iterator itr2) {
				iterator temp(itr1);
				itr1 = itr2;
				itr2 = temp;
			}
			bool operator<(iterator const& itr) const {
				return this->index_ < itr.index_;
			}
			bool operator>(iterator const& itr) const {
				return this->index_ > itr.index_;
			}
			bool operator<=(iterator const& itr) const {
				return this->index_ <= itr.index_;
			}
			bool operator>=(iterator const& itr) const {
				return this->index_ >= itr.index_;
			}
			bool operator==(iterator const& itr) const {
				return this->index_ == itr.index_;
			}
			bool operator!=(iterator const& itr) const {
				return this->index_ != itr.index_;
			}
		};
		

		iterator begin() { return iterator(this->ptr, 0); }
		iterator end() { return iterator(this->ptr, this->size); }

		class const_iterator {
		private:
			friend class Array<T>;
			using value_type = T;
			using difference_type = std::ptrdiff_t;
			using pointer = T const* ;
			using reference = T const& ;
			using iterator_category = std::random_access_iterator_tag;
			pointer ptr_;
			std::size_t index_;
		public:
			const_iterator() : ptr_(nullptr), index_(0) {}
			const_iterator(pointer ar, std::size_t index) : ptr_(ar), index_(index) {}
			const_iterator(const_iterator const& itr) : ptr_(itr.ptr_), index_(itr.index_) {}

			~const_iterator() {}

			const_iterator& operator=(const_iterator const& itr) {
				if (this == &itr)
					return *this;
				this->ptr_ = itr.ptr_;
				this->index_ = itr.index_;
				return *this;
			}

			reference operator*() {
				return this->ptr_[index_];
			}

			pointer operator->() {
				return this->ptr_;
			}

			const_iterator& operator++() {
				++this->index_;
				return *this;
			}
			const_iterator operator++(int) {
				iterator itr(this->ar, this->index_);
				this->index_++;
				return itr;
			}
			const_iterator& operator--() {
				this->index_--;
				return *this;
			}
			const_iterator operator--(int) {
				const_iterator itr(this->ar, this->index_);
				this->index_--;
				return itr;
			}
			friend const_iterator operator+(const_iterator const& itr1, const_iterator const& itr2) {
				assert(&itr1.ar == &itr2.ar);
				const_iterator itr3(itr1.ar, itr1.index_ + itr2.index_);
				return itr3;
			}
			friend const_iterator operator-(const_iterator const& itr1, const_iterator const& itr2) {
				assert(&itr1.ar == &itr2.ar);
				const_iterator itr3(itr1.ar, itr1.index_ - itr2.index_);
				return itr3;
			}
			friend const_iterator operator+(const_iterator const& itr1, difference_type num) {
				const_iterator itr2(itr1.ar, itr1.index_ + num);
				return itr2;
			}
			friend const_iterator operator+(difference_type num, const_iterator const& itr1) {
				return itr1 + num;
			}
			friend const_iterator operator-(difference_type num, const_iterator const& itr1) {
				const_iterator itr2(itr1.ptr_, itr1.index_ - num);
				return itr2;
			}
			friend const_iterator operator-(const_iterator const& itr1, difference_type num) {
				const_iterator itr2(itr1.ptr_, itr1.index_ - num);
				return itr2;
			}
			const_iterator& operator+=(difference_type num) {
				this->index_ += num;
				return *this;
			}
			const_iterator& operator-=(difference_type num) {
				this->index_ -= num;
				return *this;
			}
			reference operator[](difference_type ind) {
				return this->ptr_[ind];
			}
			friend void swap(const_iterator itr1, const_iterator itr2) {
				const_iterator temp(itr1);
				itr1 = itr2;
				itr2 = temp;
			}
			bool operator<(const_iterator const& itr) const {
				return this->index_ < itr.index_;
			}
			bool operator>(const_iterator const& itr) const {
				return this->index_ > itr.index_;
			}
			bool operator<=(const_iterator const& itr) const {
				return this->index_ <= itr.index_;
			}
			bool operator>=(const_iterator const& itr) const {
				return this->index_ >= itr.index_;
			}
			bool operator==(const_iterator const& itr) const {
				return this->index_ == itr.index_;
			}
			bool operator!=(const_iterator const& itr) const {
				return this->index_ != itr.index_;
			}
		};

		
		const_iterator begin() const { return const_iterator(this->ptr, 0); }
		const_iterator end() const { return const_iterator(this->ptr, this->size); }

		const_iterator cbegin() const { return const_iterator(this->ptr, 0); }
		const_iterator cend() const { return const_iterator(this->ptr, this->size); }
		
		// Constructors.
		Array();

		Array(std::size_t const& n);

		Array(Array const& A);

		Array(T const * const x, std::size_t const& N);

		Array(Array&& A);

		Array(Matrix<T> const& M);

		Array(Matrix<T>&& M);

		Array(std::initializer_list<T> const& L);

		// Copy assignment.
		Array& operator=(Array const& A);

		Array& operator=(Array&& A);

		Array& operator=(std::initializer_list<T> const& L);

		// Destructor.
		virtual ~Array();


		// Member functions.
		void set_size(std::size_t n);

		T& operator()(std::size_t const& n);
		T const& operator()(std::size_t const& n) const;

		Array operator-() const;

		std::size_t lenght() const { return this->size; }

		// Returns a pointer to const to the private dynamically allocated array.
		T const * get_ptr() const { return this->ptr; }

		double norm() const;

		T sum() const;

		T mean() const;

		Array<T>& add(T const& x);

		Array<T>& add(Array<T> const& A);

		Array<T>& add(std::initializer_list<T> const& L);

		void fill(T const& num);

		void reset();

		template<typename func> // Passa func come func&& o func const& per evitare copia di funzione
		void apply(func&& f,
				   typename Array<T>::iterator from = Array<T>::iterator(),
				   typename Array<T>::iterator to = Array<T>::iterator(),
				   std::size_t inc = 1);

		Array<T> each_prod(Array<T> const& A) const;
		Array<T> each_div(Array<T> const& A) const;

		// explicit type casting
		template<typename type>
		explicit operator Array<type>() const;

		Array<T> conj() const;

		friend void set_seed<T>(unsigned int seed);

		void r_uniform(T from = static_cast<T>(0.), T to = static_cast<T> (1.));

		void to_file(std::string const& fname) const;

		// Friend functions
		friend std::ostream& operator<< (std::ostream& out, Array const& A) {
			if (A.ptr == nullptr) {
				out << "[nullptr]" << std::endl;
			}
			else {
				out << "[";
				for (std::size_t i = 0; i < A.size - 1; ++i) {
					out << A(i) << ", ";
				}
				out << A(A.size - 1) << "]";
			}
			return out;
		}

		friend Array<T>(act::operator*<T>) (Matrix<T> const& M, Array<T> const& A);

	};
	template<typename T>
	std::mt19937_64 act::Array<T>::r_generator{ std::random_device{}() };
	// classic initialization std::mt19937 gen(rd());

	template<typename T>
	static void set_seed(unsigned int seed) {
		act::Array<T>::r_generator.seed(seed);
	}

	// template<typename T>
	// void act::Array<T>::set_seed(unsigned int seed) {
	// 	act::Array<T>::r_generator.seed(seed);
	// }

	template<class T> Array<T>::Array() : size(0), ptr(nullptr){}

	template<class T>
	Array<T>::Array(std::size_t const& n) : size(n), ptr(new T[n]{}) {}

	template<class T>
	Array<T>::Array(Array<T> const& A) : Array(A.size) {
		for (std::size_t i = 0; i < this->size; ++i)
			ptr[i] = A(i);
	}

	template<class T>
	Array<T>::Array(T const * const x, std::size_t const& N) : Array(N) {
		assert(x != nullptr);
		for (std::size_t i = 0; i < this->size; ++i) {
			ptr[i] = x[i];
		}
	}

	template<class T>
	Array<T>::Array(Array&& A) : size(A.size), ptr(A.ptr) {
		A.ptr = nullptr;
		A.size = 0;
	}

	template<class T>
	Array<T>::Array(Matrix<T> const& M) : size(M.rows * M.cols), ptr(new T[M.rows * M.cols]) {
		for (std::size_t i = 0; i < M.size; ++i) {
			this->ptr[i] = M.ptr[i];
		}
	}

	template<class T>
	Array<T>::Array(Matrix<T>&& M) : size(M.rows * M.cols), ptr(M.ptr) {
		M.ptr = nullptr;
		M.size = 0;
		M.cols = 0;
		M.rows = 0;
	}

	template<class T>
	Array<T>::Array(std::initializer_list<T> const& L) : Array(L.size()) {
		std::size_t c = 0;
		for (auto const& i : L) {
			ptr[c] = i;
			++c;
		}
	}

	template<class T>
	Array<T>& Array<T>::operator=(Array const& A) {
		if (this == &A) { return *this; }
		if (this->size != A.lenght()) {
			delete[] ptr;
			ptr = new T[A.lenght()];
			this->size = A.lenght();
		}
		for (std::size_t i = 0; i < this->size; ++i) {
			ptr[i] = A(i);
		}
		return *this;
	}

	template<class T>
	Array<T>& Array<T>::operator=(Array&& A) {
		this->size = A.size;
		if (this->ptr != nullptr)
			delete[] this->ptr;
		this->ptr = A.ptr;
		A.ptr = nullptr;
		A.size = 0;
		return *this;
	}

	template<class T>
	Array<T>& Array<T>::operator=(std::initializer_list<T> const& L) {
		delete[] ptr;
		ptr = new T[L.size()];
		this->size = L.size();
		std::size_t c = 0;
		for (auto const& i : L) {
			ptr[c] = i;
			++c;
		}
	}

	template<class T>
	Array<T>::~Array() { delete[] ptr; }

	template<class T>
	T& Array<T>::operator()(std::size_t const& n) { return ptr[n]; }
	template<class T>
	T const& Array<T>::operator()(std::size_t const& n) const { return ptr[n]; }

	template<class T> void Array<T>::set_size(std::size_t n) {
		assert(this->ptr == nullptr);
		this->ptr = new T[n]{};
		this->size = n;
	}

	template<class T> Array<T> Array<T>::operator-() const {
		Array<T> A(this->size);
		for (std::size_t i = 0; i < A.size; ++i) {
			A.ptr[i] = -this->ptr[i];
		}
		return A;
	}

	template<class T> double Array<T>::norm() const {
		double n{};
		for (std::size_t i = 0; i < this->size; ++i)
			n += (this->ptr[i]) * (this->ptr[i]);
		n = std::sqrt(static_cast<double>(n));
		return n;
	}

	template<class T> T Array<T>::sum() const {
		T s{};
		for (std::size_t i = 0; i < this->size; ++i) {
			s += ptr[i];
		}
		return s;
	}

	template<class T> T Array<T>::mean() const {
		return this->sum() / static_cast<double>(this->size);
	}

	template<class T> Array<T>& Array<T>::add(T const& x) {
		T *p(this->ptr);
		this->ptr = nullptr;
		this->size += 1;
		this->ptr = new T[this->size];
		for (std::size_t i = 0; i < this->size - 1; ++i)
			this->ptr[i] = p[i];
		this->ptr[this->size - 1] = x;
		delete[] p;
		p = nullptr;
		return *this;
	}

	template<class T> Array<T>& Array<T>::add(Array<T> const& A) {
		T *p(this->ptr);
		this->ptr = nullptr;
		std::size_t temp_size = this->size;
		this->size += A.size;
		this->ptr = new T[this->size];
		for (std::size_t i = 0; i < temp_size; ++i)
			this->ptr[i] = p[i];
		for (std::size_t i = 0; i < A.size; ++i)
			this->ptr[i + temp_size] = A.ptr[i];
		delete[] p;
		p = nullptr;
		return *this;
	}

	template <class T> Array<T>& Array<T>::add(std::initializer_list<T> const& L) {
		return this->add(Array<T>(L));
	}

	template<class T> void Array<T>::fill(T const& num) {
		for (std::size_t i = 0; i < this->size; ++i)
			this->ptr[i] = num;
	}

	template<class T> void Array<T>::reset() {
		this->size = 0;
		delete[] ptr;
		ptr = nullptr;
	}

	template<class T>
	template<typename func>
	void Array<T>::apply(
		func&& f, 
		typename Array<T>::iterator from,
		typename Array<T>::iterator to,
		std::size_t inc
	) {
		assert(this->ptr && "Void object.");
		if (!from.ptr_)
			from = this->begin();
		if (!to.ptr_)
			to = this->end();
		
		for (auto&& i = from; i < to; i += inc) {
			f(*i);
		}
	}

	template<class T> Array<T> Array<T>::each_prod(Array<T> const& A) const {
		assert(this->size == A.size && "Arrays have wrong sizes.");
		Array<T> res(this->size);
		for(std::size_t i = 0; i < this->size; ++i) {
			res(i) = this->ptr[i] * A.ptr[i];
		}
		return res;
	}

	template<class T> Array<T> Array<T>::each_div(Array<T> const& A) const {
		assert(this->size == A.size && "Arrays have wrong sizes.");
		Array<T> res(this->size);
		for(std::size_t i = 0; i < this->size; ++i) {
			res(i) = this->ptr[i] / A.ptr[i];
		}
		return res;
	}

	template<class T>
	template<typename type>
	Array<T>::operator Array<type>() const {
		Array<type> new_array(this->size);
		for (std::size_t i = 0; i < this->size; ++i)
			new_array(i) = static_cast<type>(this->ptr[i]);
		return new_array;
	}

	template<class T> Array<T> Array<T>::conj() const {
		Array<T> c_arr(this->size);
		for (std::size_t i = 0; i < this->size; ++i)
			c_arr(i) = std::conj(this->ptr[i]);
		return c_arr;
	}

	template<class T>
	void act::Array<T>::r_uniform(T from, T to) {
		auto&& udist = act::Uniform_distribution<T>::uniform_distribution(from, to);
		for(auto&& v : *this) {
			v = udist(act::Array<T>::r_generator);
		}
	}

	template<class T> void Array<T>::to_file(std::string const& fname) const {
		std::ofstream file(fname);
		if (this->ptr) {
			if (std::is_fundamental<T>::value ||
				std::is_same<T, std::complex<double>>::value ||
				std::is_same<T, std::complex<float>>::value) {
				file << std::setprecision(std::numeric_limits<T>::max_digits10);
				for (auto const& i : *this) {
					file << i << std::endl;
				}
			}
		}
	}

	// Save a list of Arrays into file fname.
	template<typename T>
	void to_file(std::string const& fname, std::initializer_list<act::Array<T>> L) {
		std::size_t *sizes = new std::size_t[L.size()];
		std::size_t count(0), len;
		for (auto const& el : L) {
			sizes[count] = el.lenght();
			if (count > 0)
				assert(sizes[count] == sizes[count - 1] && "Arrays have wrong dimension.");
			++count;
		}
		len = sizes[0];
		delete[] sizes;
		if (std::is_fundamental<T>::value ||
			std::is_same<T, std::complex<double>>::value ||
			std::is_same<T, std::complex<float>>::value) {
			std::ofstream file(fname);
			file << std::setprecision(std::numeric_limits<T>::max_digits10);
			for(std::size_t r = 0; r < len; ++r) {
				for(auto const& el : L) {
					file << el(r) << " ";
				}
				file << std::endl;
			}
		}
	}

	template<>
	double Array<cx<double>>::norm() const {
		double n{};
		for (std::size_t i = 0; i < this->size; ++i)
			n += std::norm(this->ptr[i]);
		n = std::sqrt(n);
		return n;
	}
	template<>
	double Array<cx<float>>::norm() const {
		float n{};
		for (std::size_t i = 0; i < this->size; ++i)
			n += std::norm(this->ptr[i]);
		n = static_cast<double>(std::sqrt(n));
		return n;
	}
	// Basic Array-Array operations
	// ============================
	template<typename S>
	Array<S> operator+(Array<S> const& A, Array<S> const& B) {
		assert(A.lenght() == B.lenght());
		Array<S> R(A.lenght());
		for (std::size_t i = 0; i < R.lenght(); ++i) {
			R(i) = A(i) + B(i);
		}
		return R;
	}

	template<typename S>
	Array<S> operator-(Array<S> const& A, Array<S> const& B) {
		assert(A.lenght() == B.lenght());
		Array<S> R(A.lenght());
		for (std::size_t i = 0; i < R.lenght(); ++i) {
			R(i) = A(i) - B(i);
		}
		return R;
	}

	template<typename S>
	S operator*(Array<S> const& A, Array<S> const& B) {
		assert(A.lenght() == B.lenght());
		S c{};
		for (std::size_t i = 0; i < A.lenght(); ++i) {
			c += A(i) * B(i);
		}
		return c;
	}

    // Array-scalar operations.
	template<typename T>
	Array<T> operator*(Array<T> const& A, T const& n) {
		Array<T> res(A);
		for(std::size_t i = 0; i < res.lenght(); ++i)
			res(i) *= n;
		return res;
	}
	template<typename T>
	Array<T> operator*(T const& n, Array<T> const& A) {
		return A*n;
	}

	template<typename T>
	Array<T> operator/(Array<T> const& A, T const& n) {
		Array<T> res(A);
		for(std::size_t i = 0; i < res.lenght(); ++i)
			res(i) /= n;
		return res;
	}
	template<typename T>
	Array<T> operator/(T const& n, Array<T> const& A) {
		return A/n;
	}


	template<typename T>
	Array<T> operator+(Array<T> const& A, T const& n) {
		Array<T> res(A);
		for(std::size_t i = 0; i < res.lenght(); ++i)
			res(i) += n;
		return res;
	}
	template<typename T>
	Array<T> operator+(T const& n, Array<T> const& A) {
		return A+n;
	}

	template<typename T>
	Array<T> operator-(Array<T> const& A, T const& n) {
		Array<T> res(A);
		for(std::size_t i = 0; i < res.lenght(); ++i)
			res(i) -= n;
		return res;
	}
	template<typename T>
	Array<T> operator-(T const& n, Array<T> const& A) {
		return A-n;
	}

	// Stack two Arrays.
	template<typename T>
	Matrix<T> stack(Array<T> const& A, Array<T> const& B) {
		assert(A.lenght() == B.lenght() && "Arrays must have same lenght.");
		Matrix<T> res(A.lenght(), 2);
		for (std::size_t r = 0; r < res.get_rows(); ++r) {
			res(r, 0) = A(r);
			res(r, 1) = B(r);
		}
		return res;
	}
	// Stack a list of Arrays.
	template<typename T>
	Matrix<T> stack(std::initializer_list<act::Array<T>> const& L) {
		std::size_t count(0);
		std::size_t *sizes = new std::size_t[L.size()];

		for (auto const& el : L) {
			sizes[count] = el.lenght();
			if (count > 0)
				assert(sizes[count] == sizes[count - 1] && "Arrays have wrong dimension.");
			++count;
		}
		Matrix<T> res(sizes[0], L.size());
		delete[] sizes;
		std::size_t c(0);
		for (auto const& el : L) {
			for (std::size_t r = 0; r < res.get_rows(); ++r) {
				res(r, c) = el(r);
			}
			++c;
		}
		return res;
	}
	// ============================


	// Template matrix class.

	template<typename S>
	Matrix<S> operator* (Matrix<S> const& M, Matrix<S> const& N);

	template<typename T>
	class Matrix : public Array<T> {
	private:
		// inherited from Array: ptr, size
		std::size_t rows, cols;
	public:
		friend class act::Array<T>;
		// Inherited public members:
		// norm(), mean(), lenght(), sum(), fill(num), get_ptr(), apply(func, from, to, inc), t_uniform(from, to), set_seed(n).

		// Cannot be used as auto i:Matrix because the functions name are not begin / end, but row_begin, row_end. For such for loops, use the simpler iterator.
		class row_iterator {
		private:
			act::Matrix<T> const *it; // Pointer to matrix.
			std::size_t index; // Current row index.
			act::Array<T> arr; // Current row contents.
		public:
			using value_type = act::Array<T>;
			using difference_type = std::ptrdiff_t;
			using pointer = act::Array<T>*;
			using reference = act::Array<T>&;
			using iterator_category = std::random_access_iterator_tag; // non mutable.

			row_iterator(){}
			row_iterator(act::Matrix<T> const& matr, std::size_t ind):
				it(&matr),
				index(ind) {
				if (this->index <= it->rows)
					arr = it->get_row(this->index);
				else
					arr.set_size(it->cols);
				}

			row_iterator(row_iterator const& itr):
				index(itr.index),
				it(itr.it),
				arr(itr.arr)
				{}

			~row_iterator(){}

			row_iterator& operator=(row_iterator const& itr) {
				if(this == &itr){
					return *this;
				}
				this->index = itr.index;
				this->it = itr.it;
				this->arr = itr.arr;
				return *this;
			}

			reference operator*() {
				return arr;
			}

			pointer operator->() {
				return &arr;
			}

			row_iterator& operator++() {
				++index;
				if (this->index <= it->rows)
					this->arr = it->get_row(index);
				else
					this->arr.set_size(it->cols);
				return *this;
			}
			row_iterator operator++(int) {
				row_iterator temp;
				temp.it = this->it;
				temp.index = this->index;
				temp.arr = this->arr;
				++index;
				return temp;
			}
			row_iterator& operator--() {
				--index;
				this->arr = it->get_row(index);
				return *this;
			}
			row_iterator operator--(int) {
				row_iterator temp;
				temp.it = this->it;
				temp.index = this->index;
				temp.arr = this->arr;
				--index;
				return temp;
			}
			friend row_iterator operator+(row_iterator const& itr1, row_iterator const& itr2) {
				assert(itr1.it == itr2.it);
				row_iterator temp(itr1.it, itr1.index + itr2.index);
				return temp;
			}
			friend row_iterator operator-(row_iterator const& itr1, row_iterator const& itr2) {
				assert(itr1.it == itr2.it);
				row_iterator temp;
				temp.index = itr1.index - itr2.index;
				temp.it = itr1.it;
				temp.arr = temp.it->get_row(temp.index);
				return temp;
			}
			friend row_iterator operator+(row_iterator const& itr1, difference_type num) {
				row_iterator temp(itr1.it, itr1.index + num);
				return temp;
			}
			friend row_iterator operator+(difference_type num, row_iterator const& itr1) {
				return itr1 + num;
			}

			friend row_iterator operator-(difference_type num, row_iterator const& itr1) {
				row_iterator temp(itr1.it, itr1.index - num);
				return temp;
			}
			friend row_iterator operator-(row_iterator const& itr1, difference_type num) {
				return num - itr1;
			}
			row_iterator& operator+=(difference_type num) {
				this->index += num;
				this->arr = this->it->get_row(this->index);
				return *this;
			}
			row_iterator& operator-=(difference_type num) {
				this->index -= num;
				this->arr = this->it->get_row(this->index);
				return *this;
			}
			reference operator[](difference_type ind) {
				this->arr = this->it->get_row(ind);
				return this->arr;
			}
			friend void swap(row_iterator itr1, row_iterator itr2) {
				row_iterator temp(itr1);
				itr1 = itr2;
				itr2 = temp;
			}
			bool operator<(row_iterator const& itr) const {
				return (this->index < itr.index);
			}
			bool operator>(row_iterator const& itr) const {
				return (this->index > itr.index);
			}
			bool operator<=(row_iterator const& itr) const {
				return (this->index <= itr.index);
			}
			bool operator>=(row_iterator const& itr) const {
				return (this->index >= itr.index);
			}
			bool operator==(row_iterator const& itr) const {
				return index == itr.index;
			}
			bool operator!=(row_iterator const& itr) const {
				return !(index == itr.index);
			}
		};

		row_iterator row_begin() const {
			return row_iterator(*this, 0);
		}
		row_iterator row_end() const {
			return row_iterator(*this, this->rows);
		}

		using iterator = typename act::Array<T>::iterator;
		using const_iterator = typename act::Array<T>::const_iterator;

		iterator begin() { return act::Array<T>::begin(); }
		iterator end() { return act::Array<T>::end(); }
		
		const_iterator begin() const { return act::Array<T>::begin(); }
		const_iterator end() const { return act::Array<T>::end(); }

		const_iterator cbegin() const { return act::Array<T>::cbegin(); }
		const_iterator cend() const { return act::Array<T>::cend(); }
		

		// Constructors.
		Matrix();

		Matrix(std::size_t const& n, std::size_t const& m);

		Matrix(Matrix<T> const& M);

		Matrix(Matrix<T>&& M);

		Matrix(act::Array<T> const& A, std::size_t r, std::size_t c);

		Matrix(act::Array<T>&& A, std::size_t r, std::size_t c);

		explicit Matrix(std::initializer_list<act::Array<T>> const& L);

		// Copy assignment
		Matrix& operator=(Matrix const& M);

		Matrix& operator=(Matrix&& M);

		virtual ~Matrix();

		void set_size(std::size_t n, std::size_t m);

		T& operator()(std::size_t const& i, std::size_t const& j);

		T const& operator()(std::size_t const& i, std::size_t const& j) const;

		Matrix operator-() const;

		T trace() const;

		Matrix adjoint() const;

		Array<T> get_diag() const;

		Array<T> get_row(std::size_t row) const;

		Array<T> get_col(std::size_t col) const;

		void reset();

		std::size_t get_rows() const;

		std::size_t get_cols() const;

		void identity();

		void identity(std::size_t n);

        // Stack the matrix mat under *this.
		Matrix<T>& vstack(act::Matrix<T>const& mat);

		// Stack the array a under *this.
		Matrix<T>& vstack(act::Array<T> const& a);

        // Stack the matrix mat to the right of *this.
		Matrix<T>& hstack(act::Matrix<T> const& mat);

		// Stack the array a to the right of *this.
		Matrix<T>& hstack(act::Array<T> const& a);

		Matrix<T> conj() const;

		Matrix<T> each_prod(Matrix<T> const& M) const;

		Matrix<T> each_div(Matrix<T> const& M) const;

		void to_file(std::string const& fname) const;

		friend std::ostream& operator<< (std::ostream& out, Matrix const& M) {
			if (M.ptr == nullptr) {
				out << "[nullptr]" << std::endl;
			}
			else {
				for (std::size_t i = 0; i < M.rows; ++i) {
					out << "[ ";
					for (std::size_t j = 0; j < M.cols; ++j) {
						out << M(i, j) << " ";
					}
					out << "]";
					out << std::endl;
				}
			}
			return out;
		}
		friend Matrix<T>(act::operator*<T>)(Matrix<T> const& M, Matrix<T> const& N);
		friend Array<T>(act::operator*<T>)(Matrix<T> const& M, Array<T> const& A);

		//explicit cast operator.
		template<typename type>
		explicit operator Matrix<type>() const;

		// Delete all the 'add' functions inherited from Array.
		Array<T>& add(Array<T> const& A) = delete;
		Array<T>& add(T const& A) = delete;
		Array<T>& add(std::initializer_list<T> const& L) = delete;
	};

	template<class T>
	Matrix<T>::Matrix() : Array<T>(), rows(0), cols(0) {}

	template<class T>
	Matrix<T>::Matrix(std::size_t const& n, std::size_t const& m): Array<T>(n * m), rows(n), cols(m) {}

	template<class T>
	Matrix<T>::Matrix(Matrix<T> const& M) : Array<T>(M.size), rows(M.rows), cols(M.cols) {
		for (std::size_t i = 0; i < this->size; ++i) {
			this->ptr[i] = M.ptr[i];
		}
	}

	template<class T>
	Matrix<T>::Matrix(Matrix<T>&& M) : Array<T>(), rows(M.rows), cols(M.cols) {
		this->size = M.size;
		this->ptr = M.ptr;
		M.ptr = nullptr;
		M.rows = 0;
		M.cols = 0;
	}

	template<class T>
	Matrix<T>::Matrix(act::Array<T> const& A, std::size_t r, std::size_t c) : Array<T>() {
		assert(r * c == A.size);
		this->rows = r;
		this->cols = c;
		this->size = r * c;
		this->ptr = new T[this->size];
		std::size_t n{};
		for (std::size_t i = 0; i < this->rows; ++i) {
			for (std::size_t j = 0; j < this->cols; ++j) {
				// this->ptr[i * this->cols + j] = A.ptr[n];
				this->operator()(i,j) = A.ptr[n];
				++n;
			}
		}
	}

	template<class T>
	Matrix<T>::Matrix(act::Array<T>&& A, std::size_t r, std::size_t c) : Array<T>() {
		assert(r * c == A.size);
		this->rows = r;
		this->cols = c;
		this->size = r * c;
		this->ptr = A.ptr;
		A.ptr = nullptr;
	}

	template<class T>
	Matrix<T>::Matrix(std::initializer_list<act::Array<T>> const& L) : Array<T>() {
		std::size_t c{ 0 }, d{ 0 };
		std::size_t * sizes = new std::size_t[L.size()];
		this->ptr = new T[L.size() * (L.begin())->size];
		this->rows = L.size();
		for (auto const& i : L) {
			sizes[c] = i.size;
			if (c > 0)
				assert(sizes[c] == sizes[c - 1] && "Arrays have wrong dimension.");
			if(c == 0)
				this->cols = i.size;
			d = 0;
			for (std::size_t j = 0; j < i.size; ++j) {
				this->ptr[c * this->cols + d] = i(j);
				++d;
			}
			++c;
		}
		this->size = this->rows * this->cols;
		delete[] sizes;
	}

	template<class T>
	Matrix<T>& Matrix<T>::operator=(Matrix const& M) {
		if (this == &M)
			return *this;
		if (this->rows != M.rows || this->cols != M.cols) {
			delete[] this->ptr;
			this->rows = M.rows;
			this->cols = M.cols;
			this->size = rows * cols;
			this->ptr = new T[this->size];
		}
		for (std::size_t i = 0; i < this->size; ++i) {
			this->ptr[i] = M.ptr[i];
		}
		return *this;
	}

	template<class T>
	Matrix<T>& Matrix<T>::operator=(Matrix&& M) {
		if (this->ptr != nullptr)
			delete[] this->ptr;
		this->ptr = M.ptr;
		this->rows = M.rows;
		this->cols = M.cols;
		this->size = M.size;
		M.ptr = nullptr;
		return *this;
	}

	template<class T>
	Matrix<T>::~Matrix() {}

	template<class T>
	void Matrix<T>::set_size(std::size_t n, std::size_t m) {
		Array<T>::set_size(n * m);
		this->rows = n;
		this->cols = m;
	}

	template<class T>
	T& Matrix<T>::operator()(std::size_t const& i, std::size_t const& j) {
		return this->ptr[i * this->cols + j];
	}

	template<class T>
	T const& Matrix<T>::operator()(std::size_t const& i, std::size_t const& j) const {
		return this->ptr[i * this->cols + j];
	}

	template<class T>
	Matrix<T> Matrix<T>::operator-() const {
		Matrix<T> M(this->rows, this->cols);
		for (std::size_t i = 0; i < M.size; ++i)
			M.ptr[i] = -this->ptr[i];
		return M;
	}

	template<class T>
	T Matrix<T>::trace() const {
		assert(this->rows == this->cols);
		T t{};
		for (std::size_t i = 0; i < this->rows; ++i) {
			t += this->ptr[i * this->cols + i];
		}
		return t;
	}

	template<class T>
	Matrix<T> Matrix<T>::adjoint() const {
		Matrix<T> M(this->rows, this->cols);
		for (std::size_t i = 0; i < this->rows; ++i) {
			for (std::size_t j = 0; j < this->cols; ++j) {
				M(i, j) = this->operator()(j, i);
			}
		}
		return M;
	}

	template<>
	Matrix<cx<double>> Matrix<cx<double>>::adjoint() const {
		Matrix<cx<double>> M(this->rows, this->cols);
		for (std::size_t i = 0; i < this->rows; ++i) {
			for (std::size_t j = 0; j < this->cols; ++j) {
				M(i, j) = std::conj(this->operator()(j, i));
			}
		}
		return M;
	}
	template<>
	Matrix<cx<float>> Matrix<cx<float>>::adjoint() const {
		Matrix<cx<float>> M(this->rows, this->cols);
		for (std::size_t i = 0; i < this->rows; ++i) {
			for (std::size_t j = 0; j < this->cols; ++j) {
				M(i, j) = std::conj(this->operator()(j, i));
			}
		}
		return M;
	}

	template<class T>
	Array<T> Matrix<T>::get_diag() const {
		assert(this->cols == this->rows);
		Array<T> R(this->rows);
		for (std::size_t i = 0; i < R.lenght(); ++i)
			R(i) = this->operator()(i, i);
		return R;
	}

	template<class T>
	Array<T> Matrix<T>::get_row(std::size_t row) const {
		Array<T> R(this->cols);
		for (std::size_t i = 0; i < this->cols; ++i) {
			R(i) = this->operator()(row, i);
		}
		return R;
	}

	template<class T>
	Array<T> Matrix<T>::get_col(std::size_t col) const {
		Array<T> R(this->rows);
		for (std::size_t i = 0; i < this->rows; ++i) {
			R(i) = this->operator()(i, col);
		}
		return R;
	}

	template<class T>
	void Matrix<T>::reset() {
		Array<T>::reset();
		this->rows = 0;
		this->cols = 0;
	}

	template<class T>
	std::size_t Matrix<T>::get_rows() const {
		return this->rows;
	}
	template<class T>
	std::size_t Matrix<T>::get_cols() const {
		return this->cols;
	}

	template<class T>
	void Matrix<T>::identity() {
		assert(this->rows == this->cols && "non-square matrix error.");
		for(std::size_t r = 0; r < this->rows; ++r) {
			for(std::size_t c = 0;  c < this->cols; ++c) {
				if(r == c)
					this->ptr[r * this->cols + c] = static_cast<T>(1.);
				else
					this->ptr[r * this->cols + c] = static_cast<T>(0.);
			}
		}
	}

	template<class T>
	void Matrix<T>::identity(std::size_t n) {
		if(this->cols != n || this->rows != n){
			delete[] this->ptr;
			this->ptr = new T[n * n]{};
			this->rows = n;
			this->cols = n;
			this->size = n * n;
		}
		for(std::size_t rr = 0; rr < this->rows; ++rr) {
			for(std::size_t cc = 0;  cc < this->cols; ++cc) {
				if(rr == cc)
					this->ptr[rr * this->cols + cc] = static_cast<T>(1.);
				else
					this->ptr[rr * this->cols + cc] = static_cast<T>(0.);
			}
		}
	}

	template<class T>
	Matrix<T>& Matrix<T>::vstack(act::Matrix<T>const& mat) {
		assert(this->cols == mat.cols && this->ptr && "unsupported size.");
		std::size_t old_rows, old_size;
		old_rows = this->rows;
		old_size = this->size;
		T* temp_ptr(this->ptr);
		this->rows = old_rows + mat.rows;
		this->size = this->rows * this->cols;
		this->ptr = new T[this->size]{};
		for(std::size_t n = 0, k = 0; n < this->size; ++n) {
			if(n < old_size) {
				this->ptr[n] = temp_ptr[n];
			}
			else {
				this->ptr[n] = mat.ptr[k];
				++k;
			}
		}
		delete[] temp_ptr;
		return *this;
	}

	template<class T>
	Matrix<T>& Matrix<T>::vstack(act::Array<T> const& a) {
		assert(this->cols == a.size && this->ptr && "unsupported size.");
		std::size_t old_rows, old_size;
		old_rows = this->rows;
		old_size = this->size;
		T* temp_ptr(this->ptr);
		this->rows += 1;
		this->size = this->cols * this->rows;
		this->ptr = new T[this->size];
		for (std::size_t n = 0; n < old_size; ++n) {
			this->ptr[n] = temp_ptr[n];
		}
		for (std::size_t n = 0; n < a.size; ++n) {
			this->ptr[n + old_size] = a.ptr[n];
		}
		delete[] temp_ptr;
		return *this;
	}

	template<class T>
	Matrix<T>& Matrix<T>::hstack(act::Matrix<T> const& mat) {
		assert(this->rows == mat.rows && this->ptr && "unsupported size.");
		std::size_t old_cols;
		old_cols = this->cols;
		T* temp_ptr(this->ptr);
		this->cols = old_cols + mat.cols;
		this->size = this->rows * this->cols;
		this->ptr = new T[this->size]{};

		for(std::size_t r = 0, old_c = 0, new_c = 0; r < this->rows; ++r) {
			old_c = 0;
			new_c = 0;
			for(std::size_t c = 0; c < this->cols; ++c) {
				if(c < old_cols) {
					this->ptr[r * this->cols + c] = temp_ptr[r *
					old_cols + old_c];
					++old_c;
				}
				else {
					this->ptr[r * this->cols + c] = mat.ptr[r *
					mat.cols + new_c];
					++new_c;
				}
			}
		}
		delete[] temp_ptr;
		return *this;
	}

	template<class T>
	Matrix<T>& Matrix<T>::hstack(act::Array<T> const& a) {
		assert(this->cols == a.size && this->ptr && "unsupported size.");
		std::size_t old_cols, old_size;
		old_cols = this->cols;
		old_size = this->size;
		T* temp_ptr(this->ptr);
		this->cols += 1;
		this->size = this->cols * this->rows;
		this->ptr = new T[this->size];

		for (std::size_t r = 0; r < this->rows; ++r) {
			for (std::size_t c = 0; c < this->cols; ++c) {
				if (c < old_cols) {
					this->ptr[r * this->cols + c] = temp_ptr[r * old_cols + c];
				}
				else {
					this->ptr[r * this->cols + c] = a.ptr[r];
				}
			}
		}
		delete[] temp_ptr;
		return *this;
	}

	template<class T>
	Matrix<T> Matrix<T>::conj() const {
		Matrix<T> c_mat(this->rows, this->cols);
		for (std::size_t i = 0; i < this->rows; ++i) {
			for (std::size_t j = 0; j < this->cols; ++j) {
				c_mat(i, j) = std::conj(this->ptr[i * this->cols + j]);
			}
		}
		return c_mat;
	}

	template<class T>
	Matrix<T> Matrix<T>::each_prod(Matrix<T> const& M) const {
		assert(this->cols == M.cols && this->rows == M.rows && "each_prod error.");
		Matrix<T> res(this->rows, this->cols);
		for(std::size_t i = 0; i < this->size; ++i) {
			res.ptr[i] = this->ptr[i] * M.ptr[i];
		}
		return res;
	}

	template<class T>
	Matrix<T> Matrix<T>::each_div(Matrix<T> const& M) const {
		assert(this->cols == M.cols && this->rows == M.rows && "each_prod error.");
		Matrix<T> res(this->rows, this->cols);
		for(std::size_t i = 0; i < this->size; ++i) {
			res.ptr[i] = this->ptr[i] / M.ptr[i];
		}
		return res;
	}

	template<class T>
	void Matrix<T>::to_file(std::string const& fname) const {
		std::ofstream file(fname);
		if (this->ptr) {
			if (std::is_fundamental<T>::value ||
				std::is_same<T, std::complex<double>>::value ||
				std::is_same<T, std::complex<float>>::value) {
				file << std::setprecision(std::numeric_limits<T>::max_digits10);
				for (int r = 0; r < this->rows; ++r) {
					for (int c = 0; c < this->cols; ++c) {
						file << this->operator()(r, c) << " ";
					}
					file << std::endl;
				}
			}
		}
	}

	template<class T>
	template<typename type>
	Matrix<T>::operator Matrix<type>() const {
		Matrix<type> new_matrix(this->rows, this->cols);
		for (std::size_t i = 0; i < this->rows; ++i) {
			for (std::size_t j = 0; j < this->cols; ++j) {
				new_matrix(i, j) = static_cast<type>(this->ptr[i * this->cols + j]);
			}
		}
		return new_matrix;
	}
	// Matrix-Matrix operations.
	template<typename T>
	Matrix<T> operator+(Matrix<T> const& A, Matrix<T> const& B) {
		assert(A.get_rows() == B.get_rows() && A.get_cols() == B.get_cols());
		Matrix<T> C(A.get_rows(), A.get_cols());
		for (std::size_t i = 0; i < C.get_rows(); ++i) {
			for (std::size_t j = 0; j < C.get_cols(); ++j) {
				C(i, j) = A(i, j) + B(i, j);
			}
		}
		return C;
	}

	template<typename T>
	Matrix<T> operator-(Matrix<T> const& A, Matrix<T> const& B) {
		assert(A.get_rows() == B.get_rows() && A.get_cols() == B.get_cols());
		Matrix<T> C(A.get_rows(), A.get_cols());
		for (std::size_t i = 0; i < C.get_rows(); ++i) {
			for (std::size_t j = 0; j < C.get_cols(); ++j) {
				C(i, j) = A(i, j) - B(i, j);
			}
		}
		return C;
	}


	template<typename T>
	Matrix<T> operator* (Matrix<T> const& M, Matrix<T> const& N) {
		assert(M.cols == N.rows && "Matrices have wrong dimensions!");
		Matrix<T> R(M.rows, N.cols);
		for (std::size_t i = 0; i < R.rows; ++i) {
			for (std::size_t j = 0; j < R.cols; ++j) {
				R(i, j) = static_cast<T>(0.);
				for (std::size_t n = 0; n < M.cols; ++n) {
					R(i, j) += M(i, n) * N(n, j);
				}
			}
		}
		return R;
	}
	template<>
	Matrix<double> operator*(Matrix<double> const& M, Matrix<double> const& N) {
		assert(M.cols == N.rows && "Matrices have wrong dimensions!");
		Matrix<double> R(M.rows, N.cols);
		double const alpha{ 1. }, beta{ 0.0 };
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M.rows, N.cols,
			M.cols, alpha, M.ptr, M.cols, N.ptr, N.cols, beta, R.ptr, R.cols);

		return R;
	}
	template<>
	Matrix<float> operator*(Matrix<float> const& M, Matrix<float> const& N) {
		assert(M.cols == N.rows && "Matrices have wrong dimensions!");
		Matrix<float> R(M.rows, N.cols);
		float const alpha{ 1.0f }, beta{ 0.0f };
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M.rows, N.cols,
			M.cols, alpha, M.ptr, M.cols, N.ptr, N.cols, beta, R.ptr, R.cols);

		return R;
	}
	template<>
	Matrix<cx<double>> operator*(Matrix<cx<double>> const& M,
		Matrix<cx<double>> const& N) {
		assert(M.cols == N.rows && "Matrices have wrong dimensions!");
		Matrix<cx<double>> R(M.rows, N.cols);
		cx<double> alpha{ 1., 0. };
		cx<double> beta{ 0., 0. };
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M.rows,
			N.cols, M.cols, &alpha, M.ptr, M.cols, N.ptr, N.cols, &beta, R.ptr, R.cols);

		return R;
	}
	template<>
	Matrix<cx<float>> operator*(Matrix<cx<float>> const& M,
		Matrix<cx<float>> const& N) {
		assert(M.cols == N.rows && "Matrices have wrong dimensions!");
		Matrix<cx<float>> R(M.rows, N.cols);
		cx<float> alpha{ 1., 0. };
		cx<float> beta{ 0., 0. };
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M.rows, N.cols,
			M.cols, &alpha, M.ptr, M.cols, N.ptr, N.cols, &beta, R.ptr, R.cols);

		return R;
	}


	// Matrix-Array operations.
	template<typename T>
	Array<T> operator*(Matrix<T> const& M, Array<T> const& A) {
		assert(M.cols == A.size);
		Array<T> R(M.rows);
		R.fill();
		for (std::size_t i = 0; i < R.size; ++i) {
			for (std::size_t j = 0; j < M.cols; ++j) {
				R(i) += M(i, j) * A(j);
			}
		}
		return R;
	}
	template<>
	Array<double> operator*(Matrix<double> const& M, Array<double> const& A) {
		assert(M.cols == A.size);
		Array<double> R(M.rows);
		double alpha{ 1. }, beta{ 0. };
		int inc{ 1 };
		cblas_dgemv(CblasRowMajor, CblasNoTrans, M.rows, M.cols, alpha,
			M.ptr, M.cols, A.ptr, inc, beta, R.ptr, inc);
		return R;
	}
	template<>
	Array<float> operator*(Matrix<float> const& M, Array<float> const& A) {
		assert(M.cols == A.size);
		Array<float> R(M.rows);
		float alpha{ 1.f }, beta{ 0.f };
		int inc{ 1 };
		cblas_sgemv(CblasRowMajor, CblasNoTrans, M.rows, M.cols, alpha,
			M.ptr, M.cols, A.ptr, inc, beta, R.ptr, inc);
		return R;
	}
	template<>
	Array<cx<double>> operator*(Matrix<cx<double>> const& M,
		Array<cx<double>> const& A) {
		assert(M.cols == A.size);
		Array<cx<double>> R(M.rows);
		cx<double> alpha{ 1., 0. }, beta{ 0., 0. };
		int inc{ 1 };
		cblas_zgemv(CblasRowMajor, CblasNoTrans, M.rows, M.cols, &alpha,
			M.ptr, M.cols, A.ptr, inc, &beta, R.ptr, inc);
		return R;
	}
	template<>
	Array<cx<float>> operator*(Matrix<cx<float>> const& M,
		Array<cx<float>> const& A) {
		assert(M.cols == A.size);
		Array<cx<float>> R(M.rows);
		cx<float> alpha{ 1.f, 0.f }, beta{ 0.f, 0.f };
		int inc{ 1 };
		cblas_cgemv(CblasRowMajor, CblasNoTrans, M.rows, M.cols, &alpha,
			M.ptr, M.cols, A.ptr, inc, &beta, R.ptr, inc);
		return R;
	}

    // Matrix-scalar operations
	template<typename T>
	Matrix<T> operator*(Matrix<T> const& M, T const& n) {
		Matrix<T> res(M.get_rows(), M.get_cols());
		for (std::size_t r = 0; r < res.get_rows(); ++r) {
			for (std::size_t c = 0; c < res.get_cols(); ++c) {
				res(r, c) = M(r, c) * n;
			}
		}
		return res;
	}

	template<typename T>
	Matrix<T> operator*(T const& n, Matrix<T> const& M) {
		return M*n;
	}

	template<typename T>
	Matrix<T> operator/(Matrix<T> const& M, T const& n) {
		Matrix<T> res(M.get_rows(), M.get_cols());
		for (std::size_t r = 0; r < res.get_rows(); ++r) {
			for (std::size_t c = 0; c < res.get_cols(); ++c) {
				res(r, c) = M(r, c) / n;
			}
		}
		return res;
	}

	template<typename T>
	Matrix<T> operator/(T const& n, Matrix<T> const& M) {
		return M/n;
	}

	template<typename T>
	Matrix<T> operator+(Matrix<T> const& M, T const& n) {
		Matrix<T> res(M.get_rows(), M.get_cols());
		for (std::size_t r = 0; r < res.get_rows(); ++r) {
			for (std::size_t c = 0; c < res.get_cols(); ++c) {
				res(r, c) = M(r, c) + n;
			}
		}
		return res;
	}

	template<typename T>
	Matrix<T> operator+(T const& n, Matrix<T> const& M) {
		return M+n;
	}

	template<typename T>
	Matrix<T> operator-(Matrix<T> const& M, T const& n) {
		Matrix<T> res(M.get_rows(), M.get_cols());
		for (std::size_t r = 0; r < res.get_rows(); ++r) {
			for (std::size_t c = 0; c < res.get_cols(); ++c) {
				res(r, c) = M(r, c) - n;
			}
		}
		return res;
	}

	template<typename T>
	Matrix<T> operator-(T const& n, Matrix<T> const& M) {
		return M-n;
	}
}
#endif
