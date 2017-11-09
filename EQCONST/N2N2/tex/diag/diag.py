from sympy import *

q1, q2, q3, q4 = symbols("q1 q2 q3 q4")
mu1, mu2, mu3 = symbols("mu1 mu2 mu3")
l1, l2 = symbols( "l1 l2" )

a00 = mu1 * l1**2 
a11 = mu2 * l2**2 * sin(q3)**2
a22 = mu2 * l2**2 
a33 = mu3

a = Matrix( [ [a00, 0, 0, 0], [0, a11, 0, 0], [0, 0, a22, 0], [0, 0, 0, a33] ])
pprint( a )

A00 = 0
A01 = -mu2 * l2**2 * cos(q2) * sin(q3) * cos(q3)
A02 = -mu2 * l2**2 * sin(q2)
A03 = 0
A10 = mu1 * l1**2
A11 = -mu2 * l2**2 * sin(q2) * sin(q3) * cos(q3)
A12 = mu2 * l2**2 * cos(q2)
A13 = 0
A20 = 0
A21 = mu2 * l2**2 * sin(q3)**2
A22 = 0
A23 = 0

A = Matrix( [ [A00, A01, A02, A03 ], [ A10, A11, A12, A13 ], [ A20, A21, A22, A23 ] ] )
pprint(A)

i00 = mu1 * l1**2 * cos(q1)**2 + mu2 * l2**2 * ( sin(q2)**2 * sin(q3)**2 + cos(q3)**2 ) + mu3 * q4**2
i01 = - mu2 * l2**2 * sin(q2) * cos(q2) * sin(q3)**2
i02 = - mu1 * l1**2 * sin(q1) * cos(q1) - mu2 * l2**2 * cos(q2) * sin(q3) * cos(q3)

i10 = i01
i11 = mu1 * l1**2 + mu2 * l2**2 * (cos(q2)**2 * sin(q3)**2 + cos(q3)**2) + mu3 * q4**2
i12 = - mu2 * l2**2 * sin(q2) * sin(q3) * cos(q3)

i20 = i02
i21 = i12
i22 = mu1 * l1**2 * sin(q1)**2 + mu2 * l2**2 * sin(q3)**2

I = Matrix( [[ i00, i01, i02 ], [ i10, i11, i12], [ i20, i21, i22 ]] )
pprint( I )

B = Matrix( [[ a00, 0, 0, 0, A00, A10, A20 ],
             [ 0, a11, 0, 0, A01, A11, A21 ],
             [ 0, 0, a22, 0, A02, A12, A22 ],
             [ 0, 0, 0, a33, A03, A13, A23 ],
             [ A00, A01, A02, A03, i00, i01, i02 ],
             [ A10, A11, A12, A13, i10, i11, i12 ],
             [ A20, A21, A22, A23, i20, i21, i22 ],
            ])
pprint(B)
print( B.eigenvals() )
