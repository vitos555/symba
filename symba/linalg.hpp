
#ifndef _SYMBA_LINALG_HPP_
#define _SYMBA_LINALG_HPP_

#include "core.hpp"

namespace symba {
namespace linalg {
using namespace core;

template<class FieldClass, std::size_t dimension, std::size_t order, class Allocators, class VectorAllocator, class DiffFormAllocator> class DiffForm;
template<class FieldClass, std::size_t dimension, class Allocators, class VectorAllocator> class Vector;
template<class FieldClass, std::size_t dimension1, std::size_t dimension2, class Allocators, class MatrixAllocator, class TransposedMatrixAllocator, class ColumnVectorAllocator, class RowVectorAllocator> class Matrix;

template<class FieldClass, std::size_t dimension, class Allocators> class VectorAllocators {
    public:
        typedef allocator<Vector<FieldClass, dimension, Allocators, VectorAllocators<FieldClass, dimension, Allocators> > > vector_allocator;
};

template<class FieldClass, std::size_t dimension, std::size_t order, class Allocators> class DiffFormAllocators {
    public:
        typedef allocator<DiffForm<FieldClass, dimension, order, Allocators, VectorAllocators<FieldClass, dimension, Allocators>, DiffFormAllocators<FieldClass, dimension, order, Allocators> > > diff_form_allocator;
        typedef allocator<DiffForm<FieldClass, dimension, order+1, Allocators, VectorAllocators<FieldClass, dimension, Allocators>, DiffFormAllocators<FieldClass, dimension, order, Allocators> > > diff_form_inc_allocator;
        typedef allocator<DiffForm<FieldClass, dimension, order-1, Allocators, VectorAllocators<FieldClass, dimension, Allocators>, DiffFormAllocators<FieldClass, dimension, order, Allocators> > > diff_form_dec_allocator;
};

template<class FieldClass, std::size_t dimension1, std::size_t dimension2, class Allocators> class MatrixAllocators {
    public:
        typedef allocator<Matrix<FieldClass, dimension1, dimension2, Allocators, MatrixAllocators<FieldClass, dimension1, dimension2, Allocators>, MatrixAllocators<FieldClass, dimension2, dimension1, Allocators>, VectorAllocators<FieldClass, dimension1, Allocators>, VectorAllocators<FieldClass, dimension2, Allocators> > > matrix_allocator;
};

}; // namespace linalg
}; // namespace symba

#endif // _SYMBA_LINALG_HPP_