(in-package :ode)

(defparameter *eps* 1.0e-14 "Used for numeric equality.")
(defparameter *step* 1.0e-4 "Used as a base stepsize.")

(defun num-equal (a b)
  (and a b (< (abs (- a b)) *eps*)))

(defun nassoc (value alist)
  (assoc value alist :test #'num-equal))

(defun swap-matrix-layout (matrix)
  (apply #'mapcar (cons #'list matrix)))

(defun central-diff-quot (f x)
  "Approximate d/dx f(x) via the central difference quotiont.
  Override *step* for different stepsize."
  (let* ((x+ (+ x *step*))
         (x- (- x *step*))
         (fx+ (funcall f x+))
         (fx- (funcall f x-))
         (formula (lambda (fx fy)  (* 0.5 (/ 1.0 *step*) (- fx fy)))))
    (if (atom fx-)
      (funcall formula fx+ fx-)
      (mapcar formula fx+ fx-))))

(defun make-central-diff-quot (f)
  "Return new function approximating the derivative of input.
  Override *step* for different stepsize."
  (lambda (x) (central-diff-quot f x)))

(defun curry-at-nth (f x n)
  "Return a new function that calls f with x modified to it's argument at the n-th
  position."
  (lambda (y)
    (let1 (tmp (modify-nth x y n))
      (funcall f tmp))))

(defun jacobian (f x)
  (let* ((dim (length x))
         (df-by-cols (loop for i below dim
                           collect (central-diff-quot
                                     (curry-at-nth f x i)
                                     (nth i x)))))
    (swap-matrix-layout df-by-cols)))

(defun make-jacobian (f)
  (lambda (x) (jacobian f x)))

(defun shorter-or-null-after (vec n)
  (or (< (length vec) n)
      (every #'zerop (nthcdr n vec))))

(defun lower-triangular-p (matrix)
  (every #'shorter-or-null-after matrix (loop for i below (length matrix) collect i)))

(defun dot-product (u v)
  (apply #'+ (mapcar #'* u v)))

(defun make-graph (f t0 t1 stepsize)
  (let ((ft0 (funcall f t0)))
    (if (> t0 t1)
      (list ft0)
      (cons ft0 (eval-f f (+ t0 stepsize) t1 stepsize)))))

