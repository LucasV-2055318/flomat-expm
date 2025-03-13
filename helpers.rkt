#lang racket

(require flomat)

; contract out
(provide    .abs!
            flomat-special-value?
            flomat-diagonal?
            exponentiate-diagonal)

(define (.abs! A)
    (define n (nrows A))
    (define m (ncols A))
    (for ([i (in-range n)])
        (for ([j (in-range m)])
            (mset! A i j (abs (ref A i j)))
        )
    )
    A
)

(define (flomat-special-value? A)
    (for/or ([i (in-range (nrows A))])
        (for/or ([j (in-range (ncols A))])
            (define value (ref A i j))
            (or (nan? value) (infinite? value))
        )
    )
)

(define (flomat-diagonal? A)
    (define n (nrows A))
    (define m (ncols A))
    (define min-dim (min n m))
    (for/and ([i (in-range min-dim)])
        (for/and ([j (in-range min-dim)])
            (when (not (= i j))
                (when (not (equal? (ref A i j) 0.0))
                    #f
                )
            )
        )
    )
)


(define (exponentiate-diagonal A)
    (define n (nrows A))
    (define m (ncols A))
    (define min-dim (min n m))
    (define result (copy-flomat A))
    (for ([i (in-range min-dim)])
        (mset! result i i (exp (ref A i i)))
    )
    result
)