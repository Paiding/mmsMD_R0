#ifndef _DSA_H
#define _DSA_H
#define DEFAULT_CAPACITY 1024

#include <stdio.h>
#include <string.h>

template <typename T> class array
{
private:
    int _size;//accurate size , n elems means _size=n
    int _capacity;
    T* _elem;
public:
    array( int c = DEFAULT_CAPACITY ){ _elem = new T[ _capacity = c ]; _size = 0;}
    array(T const * A, int lo, int hi){ copyFrom( A, lo, hi ); }
    void copyFrom(T const * A, int lo, int hi);
    //void copy(array<T> A, int lo, int hi); not safe dont use now
    void expand();
    void print();
    void put(T const & e);
    void strput(char * e);
    void remove(int rank);
    T reduce();
    int getsize();
    int getcapacity();
    void setsize(int n);
    void clear();//keep the capacity and clean all the elements
    T & operator[](int rank) const { return _elem[rank]; };
    ~array() { delete [] _elem; }
};

template <typename T> 
void array<T>::print()
{
    for (int i = 0; i < _size; i++)
    {
        printf("%d\n",_elem[i]);
    }
}

template <typename T> 
void array<T>::copyFrom(T const * A, int lo, int hi)
{//generate new array from a c 
    _elem = new T[ _capacity = 2 * ( hi - lo ) ];
    _size = 0;
    while ( lo < hi )
    {
        _elem[ _size++ ] = A[ lo++ ];
    }
    
}
/*
template <typename T> 
void array<T>::copy(array<T> A, int lo, int hi)
{//generate new array from a c 
    _elem = new T[ _capacity = 2 * ( hi - lo ) ];
    _size = 0;
    while ( lo < hi )
    {
        _elem[ _size++ ] = A[ lo++ ];
    }
    
}
*/
template <typename T> 
void array<T>::expand()
{
    if ( _size < _capacity ) return;
    _capacity = (_capacity > DEFAULT_CAPACITY) ? _capacity : DEFAULT_CAPACITY;
    T* oldElem = _elem;
    _elem = new T[ _capacity <<= 1 ];//space times 2
    for (int i = 0; i < _size; i++)
    {
        _elem[i] = oldElem[i];
    }
    delete [] oldElem;
    //expand();
}

template <typename T>
void array<T>::put(T const & e)
{
    expand();
    _elem[_size] = e;
    _size ++;
}

template <typename T>
void array<T>::strput(char * e)
{
    expand();
    int len = strlen(e);
    char* newstr = new char[len];
    for (int i = 0; i < len; i++)
    {
        newstr[i] = e[i];
    }
    _elem[_size] = newstr;
    _size ++;
}


template <typename T>
void array<T>::remove(int rank)
{
    _elem[rank] = _elem[_size - 1];
    _size --;
}

template <typename T>
T array<T>::reduce()
{
    T* sdata = new T[_size];
    for (int i = 0; i < _size; i++)
    {
        sdata[i] = _elem[i];
    }
    
    for (int i = _size; i > 1; i >>= 1)
	{
		if (i % 2 != 0)
		{
			sdata[0] += sdata[i - 1];
		}

		for (int j = 0; j < (i/2); j++)
        {   
            sdata[j] += sdata[j + i/2];
            
        }
	}
    T sum = sdata[0];
    return sum;
}

template <typename T>
int array<T>::getsize()
{
    return _size;
}

template <typename T>
int array<T>::getcapacity()
{
    return _capacity;
}

template <typename T>
void array<T>::setsize(int n)
{
    _size = n;
    expand();
}

template <typename T>
void array<T>::clear()
{
    T* oldElem = _elem;
    _elem = new T[ _capacity ];
    delete [] oldElem;
    _size = 0;
}
//template <typename T>
//T & array<T>::operator[]( int r ) const { return _elem[ r ]; }




#endif