(in-package :ode)

(defvar *eps* 1.0e-4)

(defun central-diff-quot (f)
  "Return new function approximating the derivative of input."
  (lambda (x)
    (let* ((x+ (+ x *eps*))
           (x- (- x *eps*))
           (fx+ (funcall f x+))
           (fx- (funcall f x-))
           (formula (lambda (fx fy) (* 0.5 (/ 1.0 *eps*) (- fx fy)))))
      (mapcar formula fx+ fx-))))

(defmacro make-autonomous (f)
  (let ((x (gensym)))
    `#'(lambda (,x) (conc 1 (funcall f (first ,x) (rest ,x))))))

(defvar A '((0 0) (0.5 0.5)))
(defvar b '(0.5 0.5))
(defvar c '(0 1.0))

(setq A '((0 0) (0.5 0.5)))
(setq b '(0.5 0.5))
(setq c '(0 1.0))

(setq A '((0)))
(setq b '(1))
(setq c '(0))

(defun weight-x (x k a)
  (if (not k) x (mapcar (lambda (x y z)(+ x (* y z))) x k a)))

(defun solve-rks (f x0 t0 stepsize)
  (let ((k nil))
    (let* ((scale (lambda (x) (* stepsize x)))
           (substitution-step  (lambda (a c) (mapcar scale (funcall f (+ t0 (* c stepsize)) (weight-x x0 k a))))))
      (nconc k (mapcar substitution-step A c)))))

(defun runge-kutta-single-step (f x0 t0 stepsize)
  (let* ((ks (solve-rks f x0 t0 stepsize))
         (bks (mapcar (lambda (bi k) (mapcar (lambda (ki) (* bi ki)) k)) b ks)))
    (apply #'mapcar (cons #'+ (cons x0 bks)))))


(defun runge-kutta (f x0 t0 t1 stepsize)
  (if (> t0 t1)
    (list x0)
    (let ((next (runge-kutta-single-step f x0 t0 stepsize)))
      (cons x0 (runge-kutta f next (+ t0 stepsize) t1 stepsize)))))

(defun eval-f (f t0 t1 stepsize)
  (let ((ft0 (funcall f t0)))
    (if (> t0 t1)
      (list ft0)
      (cons ft0 (eval-f f (+ t0 stepsize) t1 stepsize)))))

(defun f1 (x y)
  (list (first y) (third y) (fourth y) (second y) (* x x)))

(let* ((x0 '(1.0 1.0 1.0 1.0 1.0))
       (t0 0.0)
       (inc 0.001)
       (t1 1))
  (time (runge-kutta 'f1 x0 t0 t1 inc))
  t)
