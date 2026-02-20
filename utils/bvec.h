#ifndef __BVEC_H__
#define __BVEC_H__

#include <stdint.h>
#include <stdio.h>

#define BVEC_INIT(prefix, type)                                                                 \
                                                                                                \
    typedef struct {                                                                            \
        type *array;                                                                            \
        size_t size;                                                                            \
        size_t capacity;                                                                        \
    } prefix##_bvec_t;                                                                          \
                                                                                                \
    static inline prefix##_bvec_t *prefix##_init(size_t capacity) {                             \
        prefix##_bvec_t *bvec = (prefix##_bvec_t *)malloc(sizeof(prefix##_bvec_t));             \
        bvec->size = 0;                                                                         \
        if (capacity == 0) return bvec;                                                         \
        bvec->array = (type *)malloc(sizeof(type) * capacity);                                  \
        if (!bvec->array) { free(bvec); return NULL; }                                          \
        bvec->capacity = capacity;                                                              \
        return bvec;                                                                            \
    }                                                                                           \
                                                                                                \
    static inline void prefix##_validate(prefix##_bvec_t *bvec) {                               \
        if (bvec->size == bvec->capacity) {                                                     \
            bvec->capacity *= 2;                                                                \
            type *temp = (type *)realloc(bvec->array, sizeof(type) * bvec->capacity);           \
            if (!temp) {                                                                        \
                fprintf(stderr, "[ERROR] couldn't reallocate array\n");                         \
                exit(1);                                                                        \
            }                                                                                   \
            bvec->array = temp;                                                                 \
        }                                                                                       \
    }                                                                                           \
                                                                                                \
    static inline void prefix##_add(prefix##_bvec_t *bvec, type element) {                      \
        prefix##_validate(bvec);                                                                \
        bvec->array[bvec->size++] = element;                                                    \
    }                                                                                           \
                                                                                                \
    static inline void prefix##_clear(prefix##_bvec_t *bvec) {                                  \
        bvec->size = 0;                                                                         \
    }                                                                                           \
                                                                                                \
                                                                                                \
    static inline void prefix##_free(prefix##_bvec_t *bvec) {                                   \
        if (!bvec) return;                                                                      \
        if (bvec->array) {                                                                      \
            free(bvec->array);                                                                  \
            bvec->array = NULL;                                                                 \
        }                                                                                       \
        bvec->size = 0;                                                                         \
        free(bvec);                                                                             \
        bvec = NULL;                                                                            \
    }

#endif