;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; These are the COOT functions that are used to build a whole tree. To use 
;;; them for adding single sugars, but without auto-accept, they need to be
;;; on the outer level and that's why they are copied them from the 
;;; add-linked-residue-tree function from: 
;;; https://github.com/pemsley/coot/blob/cfd6875b62813dff5ac492c612b6c0e2fce6a8cb/scheme/gui-add-linked-cho.scm.
;;;
;;; They are somewhat modified to not do all the graphic stuff like turning and 
;;; zooming in on the protein.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; recursively add the tree
(define (process-tree parent tree proc-func imol)
    (cond
     ((null? tree) '())
     ((list? (car tree))
      (let ((part-1 (process-tree parent (car tree) proc-func imol))
	    (part-2 (process-tree parent (cdr tree) proc-func imol)))
	(cons part-1 part-2)))
     (else
      (let ((new-res (proc-func imol parent (car tree))))
	(cons new-res
	      (process-tree new-res (cdr tree) proc-func imol))))))
