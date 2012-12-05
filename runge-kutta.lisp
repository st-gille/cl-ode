(in-package :ode)

(defvar *eps* 1.0e-14)
(defvar *step* 1.0e-4)

(defun central-diff-quot (f)
  "Return new function approximating the derivative of input."
  (lambda (x)
    (let* ((x+ (+ x *step*))
           (x- (- x *step*))
           (fx+ (funcall f x+))
           (fx- (funcall f x-))
           (formula (lambda (fx fy) (* 0.5 (/ 1.0 *step*) (- fx fy)))))
      (mapcar formula fx+ fx-))))

(defun make-autonomous (f)
    (lambda (x) (list 1 (funcall f (first x) (rest x)))))

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
  (let* ((k nil)
         (scale (lambda (x) (* stepsize x)))
         (substitution-step  (lambda (a c) (mapcar scale (funcall f (+ t0 (* c stepsize)) (weight-x x0 k a))))))
    (nconc k (mapcar substitution-step A c))))

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

(time (runge-kutta-single-step (lambda (x y) '(1)) '(1.0) 0.0 0.01))

(let ((f (lambda (x y) (list (first y) (third y) (fourth y) (second y) (* x x))))
      (x0 '(1.0 1.0 1.0 1.0 1.0))
      (t0 0.0)
      (inc 0.01)
      (t1 1))
  (time (runge-kutta f x0 t0 t1 inc))
  t)

(defun num-equal (a b)
  (and a b (< (abs (- a b)) *eps*)))

(defun find-much (item alist)
  (let ((prev alist))
    (labels ((dist (other) (abs (- item (caar other))))
             (pred (list) (or (null list) (> (caar list) item)))
             (rec (cur) (if (pred (rest (setq prev cur))) prev (rec (rest cur)))))
      (rec alist))))

(defun cached-ode-solution (the-ode init-t0 init-x0 init-stepsize)
  (let ((stepsize (max *eps* init-stepsize))
        (known-values (list (list init-t0 init-x0))))
    (flet ((append-to-values (rrest vals)
             (let ((rst (cdr rrest)))
               (setf (cdr rrest) (if rst (cons vals rst) vals)))))
      (lambda (t1)
        (let* ((closest (find-much t1 known-values))
               (t0 (caar closest))
               (x0 (cadar closest)))
          (cond
            ((num-equal t0 t1)
             (print 'look-up)
             x0)
            ((<= (abs (- t0 t1)) stepsize)
             (print 'single-step)
             (append-to-values closest
               (list t1 (runge-kutta-single-step the-ode x0 t0 (- t1 t0)))))
            (t
              (print 'full-step)
              (append-to-values closest (mapcar
                                          (lambda (&rest rest) (cons (setf t0 (+ t0 stepsize)) rest))
                                          (runge-kutta the-ode x0 t0 t1 stepsize))))))))))
