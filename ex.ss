;;; Some basic functions

(define atom?
  (lambda (x)
    (and (not (pair? x)) (not (null? x)))))

(define (pos? x)
    (> x 0))

(define (neg? x)
    (not (pos? x)))

(define add1
  (lambda (x)
    (+ x 1)))

(define sub1
  (lambda (x)
    (- x 1)))

(define (average x y)
    (/ (+ x y) 2))

(define (abs x)
    (if (> 0 x)
        (- x)
        x))

(define (square x) (* x x))

(define (cube x) (* x x x))


;;; Newton's square root method 1.1.7

(define (sqrt-iter guess previous-guess x)
    (if (good-enough? guess previous-guess)
        guess
        (sqrt-iter (improve guess x) guess x)))

(define (improve guess x)
    (average guess (/ x guess)))

(define (good-enough? guess previous-guess)
    (< ( / (abs (- guess previous-guess)) guess) 0.00001))

(define (sqrt x)
    (sqrt-iter 1.0 0.0 x))

(define (cube-root-iter guess previous-guess x)
    (if (good-enough? guess previous-guess)
        guess
        (cube-root-iter (cube-root-improve guess x) guess x)))

(define (cube-root-improve guess x)
    (/ (+ (/ x (square guess)) (* 2 guess)) 3))

(define (cube-root x)
    (cube-root-iter 1.0 0.0 x))


;;; Counting change as a tree recursive procedure

(define (count-change amount)
  (cc amount 5))
(define (cc amount kinds-of-coins)
  (cond ((= amount 0) 1)
        ((or (< amount 0) (= kinds-of-coins 0)) 0)
        (else (+ (cc amount
                     (- kinds-of-coins 1))
                 (cc (- amount
                        (first-denomination kinds-of-coins))
                     kinds-of-coins)))))
(define (first-denomination kinds-of-coins)
  (cond ((= kinds-of-coins 1) 1)
        ((= kinds-of-coins 2) 5)
        ((= kinds-of-coins 3) 10)
        ((= kinds-of-coins 4) 25)
        ((= kinds-of-coins 5) 50)))


;;; Pascale's triangle

(define (pascal row col)
    (cond ((or (= col 1) (= col row)) 1)
          ((or (> col row) (< col 1) (< row 1)) 0)
          (else (+ (pascal (sub1 row) (sub1 col)) (pascal (sub1 row) col)))))


;;; Multiplying as addition

(define (remainder num div)
    0)

(define (even? n)
    (= (remainder n 2) 0))

(define (halve n)
    (/ n 2))

(define (double n)
    (* n 2))

(define (times a b)
    (cond ((= b 0) 0)
          ((and (even? a) (even? b)) (times (double a) (halve b)))
          (else (+ a (times a (- b 1))))))


;;; GCD

(define (gcd a b)
    (if (= b 0)
        a
        (gcd b (remainder a b))))

;;; Testing for primes

(define (smallest-divisor n)
    (find-divisor n 2))

(define (find-divisor n test-divisor)
    (cond ((> (square test-divisor) n) n)
          ((divides? test-divisor n) test-divisor)
          (else (find-divisor n (add1 test-divisor)))))

(define (divides? a b)
    (= (remainder b a) 0))

(define (prime? n)
    (if (= n 1) false (= n (smallest-divisor n))))


;;; Generic sum function

(define (sum term a next b)
    (if (> a b)
        0
        (+ (term a)
           (sum term (next a) next b))))

(define (sum-cubes a b)
    (sum cube a add1 b))

(define (integral f a b dx)
    (define (add-dx x) (+ x dx))
    (* (sum f (+ a (/ dx 2.0)) add-dx b)
       dx))

(define (simpsons-integral f a b n)
    (define h (/ (- b a) n))
    (define (yk k) (f (+ a (* k h))))
    (define (simpson-term k)
        (cond ((or (= k 0) (= k n))))))

(define (sum-iter term a next b)
    (define (iter a result)
        (if (> a b)
            result
            (iter (next a) (+ result (term a)))))
    (iter a 0))

(define (product term a next b)
    (if (> a b)
        1
        (* (term a)
           (product term (next a) next b))))

(define (pi-approx n)
    (* 4 ()))

(define (accumulate combiner null-value term a next b)
    (if (> a b)
        null-value
        (combiner (term a)
                  (accumulate combiner null-value term (next a) next b))))

(define (alt-sum term a next b)
    (accumulate + 0 term a next b))

(define (filtered-accumulate combiner null-value filt term a next b)
    (if (> a b)
        null-value
        (if (filt a)
            (combiner (term a)
                      (filtered-accumulate combiner null-value filt term (next a) next b))
            (combiner null-value
                      (filtered-accumulate combiner null-value filt term (next a) next b)))))

(define (sum-square-prime a b)
    (filtered-accumulate + 0 prime? square a add1 b))

            
;;; Rational numbers implementation

(define (add-rat x y)
    (make-rat (+ (* (numer x) (denom y))
                 (* (numer y) (denom x)))
              (* (denom x) (denom y))))

(define (sub-rat x y)
  (make-rat (- (* (numer x) (denom y))
               (* (numer y) (denom x)))
            (* (denom x) (denom y))))

(define (mul-rat x y)
  (make-rat (* (numer x) (numer y))
            (* (denom x) (denom y))))

(define (div-rat x y)
  (make-rat (* (numer x) (denom y))
            (* (denom x) (numer y))))

(define (equal-rat? x y)
  (= (* (numer x) (denom y))
     (* (numer y) (denom x))))

(define (normalise x)
    (let ((n (car x))
          (d (cdr x)))
        (cond ((and (neg? n) (neg? d)) (cons (- n) (- d)))
              ((neg? d) (cons (- n) (- d)))
              (else x))))

(define (make-rat n d) 
    (let ((g (gcd n d)))
        (normalise (cons (/ n g) (/ d g)))))

(define (numer x) (car x))

(define (denom x) (cdr x))

(define one-half (make-rat 1 2))

(define one-third (make-rat 1 3))

(define test1 (make-rat (- 1) (- 2)))

(define test2 (make-rat 1 (- 2)))

(define (print-rat x)
    (newline)
    (display (numer x))
    (display "/")
    (display (denom x))
    (newline))


;;; Line segments

(define (make-point x y)
    (cons x y))

(define (x-point p)
    (car p))

(define (y-point p)
    (cdr p))

(define (make-segment start-point end-point)
    (cons start-point end-point))

(define (start-segment seg)
    (car seg))

(define (end-segment seg)
    (cdr seg))

(define (midpoint-segment seg)
    (make-point (average (x-point (start-segment seg)) (x-point (end-segment seg)))
                (average (y-point (start-segment seg)) (y-point (end-segment seg)))))

(define (print-point p)
  (newline)
  (display "(")
  (display (x-point p))
  (display ",")
  (display (y-point p))
  (display ")")
  (newline))

(define test-segment (make-segment (make-point 0 0) (make-point 2 2)))

(define (length-segment seg)
    (sqrt (+ (square (- (x-point (end-segment seg)) (x-point (start-segment seg))))
             (square (- (y-point (end-segment seg)) (y-point (start-segment seg)))))))

;;; Rectangles

(define (area-rect rect)
    (* (height-rect rect) (length-rect rect)))

(define (perim-rect rect)
    (+ (double (height-rect rect)) (double (length-rect rect))))

(define (height-rect rect)
    (car (cdr rect)))

(define (length-rect rect)
    (cdr (cdr rect)))

(define (point-a-rect rect)
    (car rect))

(define (point-c-rect rect)
    (let ((p (point-a-rect rect))
          (l (length-rect rect))
          (h (height-rect rect)))
         (make-point (+ (x-point p) l) (+ (y-point p) h))))
        
(define (make-rect point length height)
    (cons point (cons length height)))

(define test-rect (make-rect (make-point 2 2) 4 2))

;;; List exercises

(define (last-pair lst)
    (list-ref lst (- (length lst) 1)))

(define (alt-last-pair lst)
    (if (eq? (cdr lst) ())
        (car lst)
        (alt-last-pair (cdr lst))))

(define (reverse lst)
    (if (null? lst)
        ()
        (append (reverse (cdr lst)) (list (car lst)))))
       
(define us-coins (list 1 25 50 10 5))
(define uk-coins (list 100 50 20 10 5 2 1 0.5))

(define (cc amount coin-values)
  (cond ((= amount 0) 1)
        ((or (< amount 0) (no-more? coin-values)) 0)
        (else
         (+ (cc amount
                (except-first-denomination coin-values))
            (cc (- amount
                   (first-denomination coin-values))
                coin-values)))))

(define (no-more? coin-values)
    (null? coin-values))

(define (first-denomination coin-values)
    (car coin-values))

(define (except-first-denomination coin-values)
    (cdr coin-values))

;;; Binary mobile

(define (make-mobile left right)
  (list left right))
 
 (define (make-branch length structure)
  (list length structure))
 
 (define (left-branch mobile)
     (car mobile))
 
 (define (right-branch mobile)
     (car (cdr mobile)))
