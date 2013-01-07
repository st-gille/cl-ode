(in-package :ode)

(define-foreign-library libnewton
                        (:unix "libnewton.so.1.0.0")
                        (t (:default "libnewton")))
(use-foreign-library libnewton)

(defcfun "lr_decomp" :int8
         (dim :uint64)
         (matrix (:pointer (:pointer :double)))
         (pivot (:pointer :uint64)))

(defcfun "lr_solve" :int8
         (dim :uint64)
         (decomp-matrix (:pointer (:pointer :double)))
         (rhs (:pointer :double))
         (solution (:pointer :double))
         (pivot (:pointer :uint64)))

(defcfun "newtons_method" :int8
         (dim :uint64)
         (x (:pointer :double))
         (func :pointer)
         (diff_func :pointer))

(defcfun "newtons_method_impicit" :int8
         (dim :uint64)
         (x (:pointer :double))
         (y (:pointer :double))
         (func :pointer)
         (diff_func :pointer))

(defcfun "print_vector" :void (dim :uint64) (x (:pointer :double)))
(defcfun "print_vector_i" :void (dim :uint64) (x (:pointer :int)))
(defcfun "print_vector_u" :void (dim :uint64) (x (:pointer :uint64)))

(defcfun "print_matrix" :void
         (rows :uint64)
         (cols :uint64)
         (x (:pointer (:pointer :double))))

(defun set-c-vector (vec vals &optional (type :double))
  "Convert list or vector to c-array"
  (loop for val being the elements of vals
        for i from 0
        do (setf (mem-aref vec type i) val)))

(defun convert-from-c-vector (size vec &optional (type :double))
  "Convert c-array to list."
  (loop for i below size
        collect (mem-aref vec type i)))

(defmethod set-c-matrix (A rows)
  "Row-wise initialization of c-style matrices."
  (loop for row being the elements of rows
        for i from 0
        do (set-c-vector (mem-aref A :pointer i) row)))

(defmethod set-c-matrix (A (rows array))
  "Row-wise initialization of c-style matrices via 2d-arrays."
  (loop for i below (array-dimension rows 0)
        do (loop for j below (array-dimension rows 1)
                 do (setf (mem-aref (mem-aref A :pointer i) :double j)
                          (aref rows i j)))))

(defmacro make-callback (name dim f &optional (convert-to-c `set-c-vector))
  `(defcallback ,name :void ((parg :pointer) (pres ::pointer))
     (let1 (arg (convert-from-c-vector ,dim parg))
           (,convert-to-c pres (funcall ,f  arg)))))

(defun newton-solver
  (nleqs len
         &optional (initial-guess (make-array len
                                              :element-type 'double-float
                                              :initial-element 0.0d0 )))
  "Solve the equations given by <nleqs> == 0. Length must be specified.
  If not specified, the <initial-guess> is taken as the zero vector of length <len>."
  (with-foreign-object (k :double len)
    (set-c-vector k initial-guess)
    (make-callback c-f
                   len
                   nleqs)
    (make-callback c-df
                   len
                   (make-jacobian nleqs)
                   set-c-matrix)
    (when (= 1 (newtons-method len k
                               (get-callback 'c-f)
                               (get-callback 'c-df)))
      (error "cant solve system"))
    (convert-from-c-vector len k)))

(defun alloc-c-matrix (A rows cols &key initial-contents (type :double))
  "Allocate enough memory to hold <rows>*<cols> objects of
  foreign type <type>. Row-wise initialization possible."
  (loop for i below (if rows rows (length initial-contents))
        do
        (setf (mem-aref A :pointer i)
              (if initial-contents
                (foreign-alloc type :initial-contents (elt initial-contents i))
                (foreign-alloc type :count cols)))))

(defun free-c-matrix (A no-of-rows)
  "Free memory allocated by alloc-c-matrix."
  (dotimes (row no-of-rows)
    (foreign-free (mem-aref A :pointer row))))

(defmacro with-matrix ((name &optional (type :double)
                             &key initial-contents
                                  (dims nil have-dims))
                       &body body)
  "Bind <name> to a foreign matrix with the specified dimensions in <body>.
  If no dimensions are specified, <initial-contents> have to be and the dimensions are inferred."
  (with-gensyms (g-dims g-contents g-rows g-cols)
    `(let* (,@(if have-dims
                `((,g-dims ,dims))
                `((,g-contents ,initial-contents)))
            (,g-rows ,(if have-dims
                        `(first ,g-dims)
                        `(length ,g-contents)))
            (,g-cols ,(if have-dims
                        `(second ,g-dims)
                        `(length (first ,g-contents)))))
       (with-foreign-object (,name :pointer ,g-rows)
         (alloc-c-matrix ,name ,g-rows ,g-cols :type ,type
                         ,@(if (not have-dims) `(:initial-contents ,g-contents)))
         ,@body
         (free-c-matrix ,name ,g-rows)))))
