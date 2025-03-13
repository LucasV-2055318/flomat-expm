#lang racket

(require flomat "helpers.rkt")

; An improved version of flomat's expm version, closely following MATLAB's expm implementation (for < 50x50 matrices)
; ATTENTION:    - This has not been tested thoroughly (yet), use at your own risk!
;               - Function assumes that matrices are < 50x50!
;               - MATLAB recomputes the diagonal when A is in Schur form, this is not implemented!
;               - There are likely some bugs...
(define (expm A)
    ; Preliminary checks
    (when (not (flomat? A))
        (error "Matrix must be a flomat"))
    (when (flomat-special-value? A)  ; custom function -> true when A contains an inf or nan value
        (error "Matrix must not contain special values"))
    (when (not (equal? (ncols A) (nrows A)))
        (error "Matrix must be square"))

    ; diagonal matrix -> exponentiate the diagonal
    (if (flomat-diagonal? A)
        (exponentiate-diagonal A)
        (exp-general A)
    )
)

(define (exp-general A)
    (define T A)

    ; Get scaling and and pade parameters
    (define-values (s m Tpowers) (scaling-and-pade-parameters T))   

    ; rescale powers of T
    (when (not (equal? s 0.0))

        (set! T (times T (/ 1 (expt 2 s))))

        (hash-set! Tpowers 2 (times (hash-ref Tpowers 2) (/ 1 (expt 2 (* 2 s)))))

        (for ([i (in-range 4 7 2)])
            (if (not (flomat-special-value? (hash-ref Tpowers i)))
                (hash-set! Tpowers i (times (hash-ref Tpowers i) (/ 1 (expt 2 (* i s)))))
                (hash-set! Tpowers i (times (hash-ref Tpowers (i - 2)) (hash-ref Tpowers 2)))
            )
        )
    )

    (define F (pade-approx T Tpowers m))
    ; squaring phase
    (for ([k (in-range s)])
        (set! F (times F F))
    )
    F
)

(define (pade-approx T Tpowers m)
    (define c (get-pade-coeffs m))
    (define n (nrows T))
    (define I (flomat-identity n))
    (define U (flomat-zeros n n))
    (define V (flomat-zeros n n))

    (if (= m 13)
        (let* ([term1 (times 
                        T
                        (plus
                            (times
                                (hash-ref Tpowers 6)
                                (plus
                                    (times (list-ref c 13) (hash-ref Tpowers 6))
                                    (times (list-ref c 11) (hash-ref Tpowers 4))
                                    (times (list-ref c 9) (hash-ref Tpowers 2))
                                )
                            )
                            (times (list-ref c 7) (hash-ref Tpowers 6))
                            (times (list-ref c 5) (hash-ref Tpowers 4))
                            (times (list-ref c 3) (hash-ref Tpowers 2))
                            (times (list-ref c 1) I)
                        )
                    )]
                [term2  (plus
                            (times
                                (hash-ref Tpowers 6)
                                (plus
                                    (times (list-ref c 12) (hash-ref Tpowers 6))
                                    (times (list-ref c 10) (hash-ref Tpowers 4))
                                    (times (list-ref c 8) (hash-ref Tpowers 2))
                                )
                            )
                            (times (list-ref c 6) (hash-ref Tpowers 6))
                            (times (list-ref c 4) (hash-ref Tpowers 4))
                            (times (list-ref c 2) (hash-ref Tpowers 2))
                            (times (list-ref c 0) I)
                        )])
            (set! U term1)
            (set! V term2))

        ; MATLAB: numel(Tpowers) = 6 (+ 2), but this is static, so can be replaced with 6
        (let* ([strt 8])
            
            ;; Generate the powers of T in Tpowers hash
            (for ([k (in-range strt m 2)])
                (hash-set! Tpowers k (times (hash-ref Tpowers (- k 2)) (hash-ref Tpowers 2))))

            ;; Initialize U and V (Adjust for 0-based indexing by subtracting 1 from the index)
            (set! U (times (list-ref c 1) I))  ;; c(2) in 1-based becomes c(1) in 0-based
            (set! V (times (list-ref c 0) I))  ;; c(1) in 1-based becomes c(0) in 0-based

            ;; Compute U and V using the loop over coefficients (adjusting for 0-based)
            (for ([j (in-range m 1 -2)])
                (set! U (plus U (times (list-ref c j) (hash-ref Tpowers (- j 1)))))
                (set! V (plus V (times (list-ref c (- j 1)) (hash-ref Tpowers (- j 1))))))
            ;; Multiply U by T
            (set! U (times T U))
            ;;  Return U and V (if you need them)
        )
    )
    (plus (mldivide (minus V U) (times 2 U)) I)
)


(define (ell T coeff m_val)
    ; Scaling the matrix T by the coefficient
    (define scaledT (times (expt coeff (/ 1 (+ (* 2 m_val) 1))) (.abs! (copy-flomat T))))
    
    ; Compute alpha
    ; In MATLAB: normAm just uses 1-norm explicitly if dim < 50
    ; For Matrix-Herbie I assume we will not use that large matrices
    (define alpha (/ (norm (power scaledT (+ (* 2 m_val) 1)) 1) (norm T 1)))

    ; ; Compute t
    (define t (max (ceiling (/ (log (* 2 (/ alpha 2.2204e-16)) 2) (* 2 m_val))) 0))
    
    t)

(define (scaling-and-pade-parameters T)
    (define coeff
        (list (/ 1 100800)
            (/ 1 10059033600)
            (/ 1 4487938430976000)
            (/ 1 5914384781877411840000)
            (/ 1 113250775606021113483283660800000000)))

    (define theta
        (list 
            ; 3.650024139523051e-008 ; Commented out
            ; 5.317232856892575e-004 ; Commented out
            1.495585217958292e-002  ; m_vals = 3
            ; 8.536352760102745e-002 ; Commented out
            2.539398330063230e-001  ; m_vals = 5
            ; 5.414660951208968e-001 ; Commented out
            9.504178996162932e-001  ; m_vals = 7
            ; 1.473163964234804e+000 ; Commented out
            2.097847961257068e+000  ; m_vals = 9
            ; 2.811644121620263e+000 ; Commented out
            ; 3.602330066265032e+000 ; Commented out
            ; 4.458935413036850e+000 ; Commented out
            5.371920351148152e+000)) ; m_vals = 13

    (define Tpowers (make-hash))
    (hash-set! Tpowers 2 (times T T))
    (hash-set! Tpowers 4 (times (hash-ref Tpowers 2) (hash-ref Tpowers 2)))
    (hash-set! Tpowers 6 (times (hash-ref Tpowers 2) (hash-ref Tpowers 4)))
    
    (define d4 (expt (norm (hash-ref Tpowers 4) 1) 1/4))
    (define d6 (expt (norm (hash-ref Tpowers 6) 1) 1/6))

    (define eta1 (max d4 d6))
    (define s 0)

    (if (and (<= eta1 (list-ref theta 0)) (equal? (ell T (list-ref coeff 0) 3) 0.0))
        (values s 3 Tpowers)
        (if (and (<= eta1 (list-ref theta 1)) (equal? (ell T (list-ref coeff 1) 5) 0.0))
            (values s 5 Tpowers)
            (let* (  [d8 (expt (norm (times (hash-ref Tpowers 4) (hash-ref Tpowers 4)) 1) 1/8)]
                    [eta3 (max d6 d8)])
                ; MATLAB: for "small matrices"... again I assume matrices are < 50x50
                (if (and (<= eta3 (list-ref theta 2)) (equal? (ell T (list-ref coeff 2) 7) 0.0))
                    (values s 7 Tpowers)
                    (if (and (<= eta3 (list-ref theta 3)) (equal? (ell T (list-ref coeff 3) 9) 0.0))
                        (values s 9 Tpowers)
                        (let* ( [d10  (expt (norm (times (hash-ref Tpowers 4) (hash-ref Tpowers 6)) 1) 1/10)]
                                [eta4 (max d8 d10)]
                                [eta5 (min eta3 eta4)])
                            ; again for small matrices...

                            (set! s (max (ceiling (log (/ eta5 (list-ref theta 4)) 2)) 0.0))
                            (set! s (+ s (ell (times T (/ 1 (expt 2 s))) (list-ref coeff 4) 13)))

                            (when (infinite? s)
                                (set! s (log (/ (norm T 1) (list-ref theta 4)) 2))
                            )

                            (values s 13 Tpowers)
                        )
                    )
                )
            )
        )
    )
)

(define (get-pade-coeffs m)
    (cond
    [(= m 3) '(120 60 12 1)]
    [(= m 5) '(30240 15120 3360 420 30 1)]
    [(= m 7) '(17297280 8648640 1995840 277200 25200 1512 56 1)]
    [(= m 9) '(17643225600 8821612800 2075673600 302702400 30270240
                2162160 110880 3960 90 1)]
    [(= m 13) '(64764752532480000 32382376266240000 7771770303897600
                1187353796428800 129060195264000 10559470521600
                670442572800 33522128640 1323241920 40840800 960960
                16380 182 1)])
)

(define A (matrix '((0.017805101599905476 0.1722176715660912) (-0.2029362425481171 0.06295344181270353))))
(displayln (expm A))