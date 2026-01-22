#ifndef NANOFLANN_POOL_ALLOCATOR_H
#define NANOFLANN_POOL_ALLOCATOR_H

namespace nanoflann {
    /** @addtogroup memalloc_grp Memory allocation
     * @{ */

     /**
      * Pooled storage allocator
      *
      * The following routines allow for the efficient allocation of storage in
      * small chunks from a specified pool.  Rather than allowing each structure
      * to be freed Solutionly, an entire pool of storage is freed at once.
      * This method has two advantages over just using malloc() and free().  First,
      * it is far more efficient for allocating small objects, as there is
      * no overhead for remembering all the information needed to free each
      * object or consolidating fragmented memory.  Second, the decision about
      * how long to keep an object is made at the time of allocation, and there
      * is no need to track down all the objects to free them.
      *
      */

    template <typename T>
    class PooledAllocator
    {
        std::list<std::unique_ptr<T>> data;

        void internal_init()
        {
            data.clear();
        }

    public:

        /**
            Default constructor. Initializes a new pool.
         */
        PooledAllocator() { internal_init(); }

        /**
         * Destructor. Frees all the memory allocated in this pool.
         */
        ~PooledAllocator() { free_all(); }

        /** Frees all allocated memory chunks */
        void free_all()
        {
            internal_init();
        }

        T *allocate()
        {
            data.emplace_back(new T);
            return data.back().get();
        }
    };
    /** @} */
}

#endif // !NANOFLANN_POOL_ALLOCATOR_H
