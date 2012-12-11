(in-package :ode)

(defmacro aif (test-form then-form &optional else-form)
  "Anaphoric if. See Let over Lambda."
  `(let ((it ,test-form))
        (if it ,then-form ,else-form)))

(defmacro let1 (binding &body body)
  "Save two parens when using only one binding. From Land of Lisp."
  `(let (,binding)
     ,@body))

(defun nest (lst)
  "Convert a list to a nested list, like
  (nest '(a b c d)) ==> (A (B (C D)))."
  (if (null (rest lst))
    (car lst)
    (list (car lst) (nest (rest lst)))))

(defun apply-predicate (object pred)
  (if (atom pred)
    (list pred object)
    (nest (append pred (list object)))))

(defmacro define-for-predicates (name cnd doc)
  `(defmacro ,name (object &rest predicates)
     ,doc
     (once-only (object)
       (let1 (tests (mapcar (curry #'apply-predicate object) predicates))
         `(,',cnd ,@tests)))))

(define-for-predicates all-p and
  "Check if every predicate in <predicates> is true for <object>")

(define-for-predicates any-p or
  "Check if any predicate in <predicates> is true for <object>")

(defmacro with-gensym-if-not (cond (&rest bindings) &body body)
  (if cond
    `(let ,(loop for (var init-form) in bindings collect `(,var (,init-form))) ,@body)
    `(with-gensyms ,(mapcar #'car bindings) ,@body)))

(defmacro wrap-binds (form bindings &body body)
  "Nest all bindings in <bindings> within <form>."
  (if bindings
    `(,form ,(first bindings)
            (wrap-binds ,form ,(rest bindings) ,@body))
    `(progn ,@body)))

(defmacro make-alist-calling (fname &rest rest)
  "Construct an alist of symbols and arguments,
  calling the function <fname> on each argument."
  `(list ,@(mapcar (lambda (x)
                     `(cons (quote ,(car x)) (,fname ,@(second x))))
                   rest)))

(defun modify-nth (lst new-val n)
  "Make a list with a fresh cell at the n-th position."
  ;This implies that all previous cells are fresh, too.
  (cond
    ((null lst) (error "List too short."))
    ((< n 0) (error "Cannot access negative positions."))
    ((= n 0) (cons new-val (cdr lst)))
    (t (cons (car lst) (modify-nth (cdr lst) new-val (1- n))))))

(defun make-simple-list (dim &optional (formula #'+))
  "Make a list of length <dim> with elements of type double-float determined by <formula>.
  The position in the list is passed to <formula>."
  (loop for j below dim collect (coerce (funcall formula j) 'double-float)))

(defun make-simple-matrix (rows cols &optional (formula #'+))
  "Make a list of length <rows>, each element being a list of length <cols>,
  which in turn has elements of type double-float determined by <formula>.
  The current row and column is passed to <formula>."
  (loop for i below rows collect (loop for j below cols collect (coerce (funcall formula i j) 'double-float))))

(defun list-to-2d-array (list)
  "Due to http://stackoverflow.com/questions/9549568/."
  (make-array (list (length list)
                    (length (first list)))
              :initial-contents list))

(defun take-n (n lst)
  (loop for i below n
        for item in lst
        collect item))

(defun split-n (n lst)
  "Split in parts of length <n>."
  (if lst
    (cons (take-n n lst)
          (split-n n (nthcdr n lst)))))

(defun insert-at (place-to-put new-points)
  (let1 (tail (cdr place-to-put))
    (setf (cdr place-to-put) (append new-points tail))))

(defparameter *eps* 1.0e-14 "Used for numeric equality.")
(defparameter *stepsize* 1.0e-4 "Used as a base stepsize.")

(defun num-equal (a b &optional (eps *eps*))
  (and a b (< (abs (- a b)) eps)))

(defun nassoc (value alist)
  (assoc value alist :test #'num-equal))

(defun swap-matrix-layout (matrix)
  (apply #'mapcar (cons #'list matrix)))

(defun central-diff-quot (f x &optional (stepsize *stepsize*))
  "Approximate d/dx f(x) via the central difference quotiont."
  (let* ((x+ (+ x stepsize))
         (x- (- x stepsize))
         (fx+ (funcall f x+))
         (fx- (funcall f x-))
         (formula (lambda (fx fy)  (* 0.5 (/ 1.0 stepsize) (- fx fy)))))
    (if (atom fx-)
      (funcall formula fx+ fx-)
      (mapcar formula fx+ fx-))))

(defun make-central-diff-quot (f &optional (stepsize *stepsize*))
  "Return new function approximating the derivative of input."
  (lambda (x) (central-diff-quot f x stepsize)))

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

(defun make-graph (f t0 t1 &optional (stepsize *stepsize*))
  (let ((ft0 (funcall f t0)))
    (if (> t0 t1)
      (list ft0)
      (cons ft0 (make-graph f (+ t0 stepsize) t1 stepsize)))))

