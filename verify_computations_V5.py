def no_interaction(M, num):
 "Tells if the matrix satisfies the condtions of Tartar, in some sense: no interaction btw particules with same speed"
 cpt = 0
# indices = [(0, 4), (4, 0), (0, 5), (5, 0), (4, 5), (5, 4), (6, 10), (10, 6), (6, 11), (11, 6), (10, 11), (11, 10)]

 indices = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 7), (8, 8), (9, 9), (10, 10), (11, 11), (0, 4), (4, 0), (0, 6), (6, 0), (4, 5), (5, 4), (6, 10), (10, 6), (6, 11), (11, 6), (10, 11), (11, 10)]

 for ind in indices:
     #if simplify(M[ind]) != 0:
     if M[ind] !=0:
         cpt += 1
         print(num)
         print(ind)
 
 return(cpt == 0);

######################################

def hat(v):
 "Returns the corresponding 3 x 3 skew-symetric matrix 'v hat'"
 return(Matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]));

######################################

from sympy import *
init_printing(use_unicode=True)

# parameters
E = symbols("E", postive = True) # Young's modulus
G = symbols("G", postive = True) # sheat modulus
a = symbols("a", postive = True) # cross section area
k2 = symbols("k_2", postive = True); k1 = symbols("k_1", postive = True) # correction factors
k3 = symbols("k_3", postive = True)
I2 = symbols("I2", postive = True); I3 = symbols("I3", postive = True) # area moments of inertia
rho = symbols("rho", postive = True) # density in reference config.
J = Matrix([[k1*(I2+I3), 0, 0], [0, I2, 0], [0, 0, I3]])

# A
Q = diag(E/(sqrt(rho)**2), k2*G/(sqrt(rho)**2), k3*G/(sqrt(rho)**2), G/(sqrt(rho)**2), E/(sqrt(rho)**2), E/(sqrt(rho)**2))
A = Matrix(BlockMatrix([[zeros(6, 6), -Q], [-eye(6), zeros(6, 6)]]))
#print(A.eigenvals())

# B tilda
Upsilon01, Upsilon02, Upsilon03 = symbols("Upsilon_01 Upsilon_02 Upsilon_03")
Upsilon0 = Matrix([[Upsilon01], [Upsilon02], [Upsilon03]])
Upsilon0Hat = hat(Upsilon0)
Gamma01, Gamma02, Gamma03 = symbols("Gamma_01 Gamma_02 Gamma_03")
Gamma0 = Matrix([[Gamma01], [Gamma02], [Gamma03]])
Gamma0Hat = hat(Gamma0)
e1Hat = hat(Matrix([[1], [0], [0]]))

lam7, lam8, lam9, lam10 = symbols("lambda_7 lambda_8, lambda_9 lambda_10", positive = True)

M1 = rho*a*diag(lam7**2, lam8**2, lam9**2)
M2 = rho*J*diag(lam10**2, lam7**2, lam7**2)

Btild11 = zeros(6, 6)
Btild21 = Matrix(BlockMatrix([[Upsilon0Hat, Gamma0Hat + e1Hat], [zeros(3, 3), Upsilon0Hat]]))
Bta = 1/(rho*a)*Upsilon0Hat*M1
Btc = 1/rho*J.inv()*(Gamma0Hat + e1Hat)*M1
Btd = 1/rho*J.inv()*Upsilon0Hat*M2
Btild12 = Matrix(BlockMatrix([[Bta, zeros(3, 3)], [Btc, Btd]]))
Btild22 = zeros(6, 6)
Btild = -Matrix(BlockMatrix([[Btild11, Btild12], [Btild21, Btild22]]))
#print(Btild.eigenvals())

Btild0 = Btild.applyfunc(lambda x : x.subs([(Upsilon01, 0), (Upsilon02, 0), (Upsilon03, 0), (Gamma01, 0), (Gamma02, 0), (Gamma03, 0)]))
normBtild0 = Btild0.norm(2)


# L and R
Dpos = diag(lam7, lam8, lam9, lam10, lam7, lam7)
Dneg = -Dpos
invDpos = diag(1/lam7, 1/lam8, 1/lam9, 1/lam10, 1/lam7, 1/lam7)
invDneg = -invDpos
# Contrary to int the notes, here L's COLUMNS are the left-eigenvect. and R is the inverse of L^T :
L = Matrix(BlockMatrix([[eye(6), eye(6)], [Dpos, Dneg]]))
R = 1/2*Matrix(BlockMatrix([[eye(6), eye(6)], [invDpos, invDneg]]))
D = Matrix(BlockMatrix([[Dneg, zeros(6, 6)], [zeros(6, 6), Dpos]]))


# G^k tilda
Gtild1 = zeros(12, 12)
Gtild1[1, 5] = Gtild1[5, 1] = 1/2
Gtild1[2, 4] = Gtild1[4, 2] = -1/2
Gtild1[8, 10] = Gtild1[10, 8] = lam9**2/2
Gtild1[7, 11] = Gtild1[11, 7] = -lam8**2/2

Gtild2 = zeros(12, 12)
Gtild2[2, 3] = Gtild2[3, 2] = 1/2
Gtild2[0, 5] = Gtild2[5, 0] = -1/2
Gtild2[6, 11] = Gtild2[11, 6] = lam7**2/2
Gtild2[9, 8] = Gtild2[8, 9] = - lam9**2/2

Gtild3 = zeros(12, 12)
Gtild3[0, 4] = Gtild3[4, 0] = 1/2
Gtild3[1, 3] = Gtild3[3, 1] = -1/2
Gtild3[7, 9] = Gtild3[9, 7] = lam8**2/2
Gtild3[6, 10] = Gtild3[10, 6] = -lam7**2/2

Gtild4 = zeros(12, 12)
Gtild4[4, 5] = Gtild4[5, 4] = (I2 - I3)/(k1*(I2 + I3)*2)
Gtild4[10, 11] = Gtild4[11, 10] = lam7**2*(I3 - I2)/(k1*(I2 + I3)*2)
Gtild4[7, 8] = Gtild4[8, 7] = (lam9**2 - lam8**2)*a/(k1*(I2 + I3)*2)

Gtild5 = zeros(12, 12)
Gtild5[3, 5] = Gtild5[5, 3] = (I3 - k1*(I2 + I3))/(I2*2)
Gtild5[9, 11] = Gtild5[11, 9] = (k1*(I2 + I3)*lam10**2 - (lam7**2)*I3)/(I2*2)
Gtild5[6, 8] = Gtild5[8, 6] = (lam7**2 - lam9**2)*a/(I2*2)

Gtild6 = zeros(12, 12)
Gtild6[3, 4] = Gtild6[4, 3] = (k1*(I2 + I3) - I2)/(I3*2) 
Gtild6[9, 10] = Gtild6[10, 9] = ((lam7**2)*I2 - (lam10**2)*k1*(I2 + I3))/(I3*2)
Gtild6[6, 7] = Gtild6[7, 6] = (lam8**2 - lam7**2)*a/(I3*2)

Gtild7 = zeros(12, 12)
Gtild7[2, 10] = Gtild7[10, 2] = 1/2
Gtild7[1, 11] = Gtild7[11, 1] = -1/2
Gtild7[4, 8] = Gtild7[8, 4] = -1/2
Gtild7[5, 7] = Gtild7[7, 5] = 1/2

Gtild8 = zeros(12, 12)
Gtild8[0, 11] = Gtild8[11, 0] = 1/2
Gtild8[2, 9] = Gtild8[9, 2] = -1/2
Gtild8[5, 6] = Gtild8[6, 5] = -1/2
Gtild8[3, 8] = Gtild8[8, 3] = 1/2

Gtild9 = zeros(12, 12)
Gtild9[1, 9] = Gtild9[9, 1] = 1/2
Gtild9[0, 10] = Gtild9[10, 0] = -1/2
Gtild9[3, 7] = Gtild9[7, 3] = -1/2
Gtild9[4, 6] = Gtild9[6, 4] = 1/2

Gtild10 = zeros(12, 12)
Gtild10[5, 10] = Gtild10[10, 5] = 1/2
Gtild10[4, 11] = Gtild10[11, 4] = -1/2

Gtild11 = zeros(12, 12)
Gtild11[3, 11] = Gtild11[11, 3] = 1/2
Gtild11[5, 9] = Gtild11[9, 5] = -1/2

Gtild12 = zeros(12, 12)
Gtild12[4, 9] = Gtild12[9, 4] = 1/2
Gtild12[3, 10] = Gtild12[10, 3] = -1/2

# B
B = Transpose(L)*Btild*R
#print(B.eigenvals())

# B with initialy straight
B0 = B.applyfunc(lambda x : x.subs([(Upsilon01, 0), (Upsilon02, 0), (Upsilon03, 0), (Gamma01, 0), (Gamma02, 0), (Gamma03, 0)]))
#normB0 = B0.norm(2)

#q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12 = symbols("q_1, q_2, q_3, q_4, q_5, q_6, q_7, q_8, q_9, q_10, q_11, q_12")
#Q_bis = diag(q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12)
#Theta_tild_0 = Transpose(B0)*Q_bis + Q_bis*B0

# G^k
# Here, if k \leq 6 then Gk = R^T (Gtildk + lambdak Gtild(k+6)) R
G1 = (Transpose(R)*(Gtild1 + lam7*Gtild7)*R).applyfunc(simplify)
G2 = (Transpose(R)*(Gtild2 + lam8*Gtild8)*R).applyfunc(simplify)
G3 = (Transpose(R)*(Gtild3 + lam9*Gtild9)*R).applyfunc(simplify)
G4 = (Transpose(R)*(Gtild4 + lam10*Gtild10)*R).applyfunc(simplify)
G5 = (Transpose(R)*(Gtild5 + lam7*Gtild11)*R).applyfunc(simplify)
G6 = (Transpose(R)*(Gtild6 + lam7*Gtild12)*R).applyfunc(simplify)

# Here, if k > 6 then Gk = R^T (Gtild(k-6) + lambdak Gtildk) R
G7 = (Transpose(R)*(Gtild1 - lam7*Gtild7)*R).applyfunc(simplify)
G8 = (Transpose(R)*(Gtild2 - lam8*Gtild8)*R).applyfunc(simplify)
G9 = (Transpose(R)*(Gtild3 - lam9*Gtild9)*R).applyfunc(simplify)
G10 = (Transpose(R)*(Gtild4 - lam10*Gtild10)*R).applyfunc(simplify)
G11 = (Transpose(R)*(Gtild5 - lam7*Gtild11)*R).applyfunc(simplify)
G12 = (Transpose(R)*(Gtild6 - lam7*Gtild12)*R).applyfunc(simplify)

# To check Tartar's condition:
#print(no_interaction(G1, 1), no_interaction(G2, 2), no_interaction(G3, 3), no_interaction(G4, 4), no_interaction(G5, 5), no_interaction(G6, 6), no_interaction(G7, 7), no_interaction(G8, 8), no_interaction(G9, 9), no_interaction(G10, 10), no_interaction(G11, 11), no_interaction(G12, 12))

# To export the matrices G^k:
all_mat_tex = ""
for k in range(1, 13):
 mat = eval("((G%s*8).applyfunc(lambda x : collect(x, [lam7*lam10]))).applyfunc(simplify)" % str(k))
 mat1_str = latex(mat[:, 0:6])
 mat2_str = latex(mat[:, 6:12])
 mat1_str = mat1_str.replace("\end{matrix}\\right]", "\end{matrix}\\right.")
 mat2_str = mat2_str.replace("\\left[\\begin{matrix}", "\\left.\\begin{matrix}")
 mat_tex = "$G_{%s}$ computed by sympy: \n \\begin{align*} \n \\frac{1}{8} & %s \\\ \n & %s \n \end{align*} \n  \n" % (str(k),mat1_str, mat2_str)
 mat_tex = mat_tex.replace("1.0", " 1 "); mat_tex = mat_tex.replace("2.0", " 2 ")
 mat_tex = mat_tex.replace("{ 1  a ", "{  a ")
 mat_tex = mat_tex.replace("\frac{ 1  \lambda", "\frac{  \lambda")
 mat_tex = mat_tex.replace("\\left(I_{2} + I_{3} - k_{1} \\left(I_{2} + I_{3}\\right)\\right)\\right)", "\\left(1 - k_{1}\\right) \\left(I_{2} + I_{3}\\right)\\right)")
 #mat_tex = mat_tex.replace("-   -", " + "), mat_tex = mat_tex.replace("", "")
 all_mat_tex += mat_tex

with open("OutputG.txt", "w") as text_file:
 text_file.write(all_mat_tex)


# To export the matrices \widetilde{B} and B
Btex = ""
tex_Btild1 = latex(Btild[:, 0:6])
tex_Btild2 = latex(Btild[:, 6:12])
tex_Btild1 = tex_Btild1.replace("\end{matrix}\\right]", "\end{matrix}\\right.")
tex_Btild2 = tex_Btild2.replace("\\left[\\begin{matrix}", "\\left.\\begin{matrix}")
Btex += "$\widetilde{B}$ computed by sympy: \n \\begin{align*} \n  & %s \n & %s \n \end{align*} \n  \n" % (tex_Btild1, tex_Btild2)

tex_B1 = latex(B[:, 0:6])
tex_B2 = latex(B[:, 6:12])
tex_B1 = tex_B1.replace("0.5", " ")
tex_B2 = tex_B2.replace("0.5", " ")
tex_B1 = tex_B1.replace("\end{matrix}\\right]", "\end{matrix}\\right.")
tex_B2 = tex_B2.replace("\\left[\\begin{matrix}", "\\left.\\begin{matrix}")
Btex += "\n \n $B$ computed by sympy: \n \\begin{align*} \n  \\frac{1}{2} & %s \\\ \n & %s \n \end{align*} \n  \n $\widetilde{B}$ for initially straight beam: \n \n \\begin{align*} \n %s \n \end{align*} \n \n $B$ for initially straight beam: \n \n  \\begin{align*} \n %s \n \end{align*}" % (tex_B1, tex_B2, latex(Btild0), latex(B0))


with open("OutputB.txt", "w") as text_file:
 text_file.write(Btex)


# To export the developped right hand side:
y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 = symbols("y_1 y_2 y_3 y_4 y_5 y_6 y_7 y_8 y_9 y_10 y_11 y_12")
y = Matrix([[y1], [y2], [y3], [y4], [y5], [y6], [y7], [y8], [y9], [y10], [y11], [y12]])
Ry = zeros(12, 1)

v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12 = symbols("v_1 v_2 v_3 v_4 v_5 v_6 v_7 v_8 v_9 v_10 v_11 v_12")
v = Matrix([[v1], [v2], [v3], [v4], [v5], [v6], [v7], [v8], [v9], [v10], [v11], [v12]])
Rv = zeros(12, 1)

Ry0 = zeros(12, 1); Rv0 = zeros(12, 1)
for k in range (12):
 Ry[k, 0] = -Btild[k, :]*y + Transpose(y)*(eval("Gtild%s" % str(k+1))*y)
 Rv[k, 0] = simplify(-B[k, :]*v + Transpose(v)*(eval("G%s" % str(k+1))*v))
 Ry0[k, 0] = -Btild0[k, :]*y + Transpose(y)*(eval("Gtild%s" % str(k+1))*y)
 Rv0[k, 0] = simplify(-B0[k, :]*v + Transpose(v)*(eval("G%s" % str(k+1))*v))

tex_right = ""
tex_right += "Developped $-\widetilde{B}y +  \widetilde{g}(y):$ \n \n \\begin{align*} \n %s \n \end{align*}" % latex(Ry)
tex_right += "\newpage Developped $-By + g(y):$\n \n \\begin{align*} \n %s \n \end{align*}" % latex(Rv)
tex_right += "\newpage Developped $-\widetilde{B}y +  \widetilde{g}(y)$ (initially straight beam): \n \n \\begin{align*} \n %s \n \end{align*}" % latex(Ry0)
tex_right += "\newpage Developped $-By + g(y)$ (initially straight beam): \n \n \\begin{align*} \n %s \n \end{align*}" % latex(Rv0)

with open("Output_right.txt", "w") as text_file:
 text_file.write(tex_right)


# To export the matrices \widetilde{G}^k:
 all_mat_tex = ""
for k in range(1, 13):
 mat = eval("Gtild%s*2" % str(k))
 mat_str = latex(mat)
 mat_str = mat_str.replace("k_{1} \\left(I_{2} + I_{3}\\right)", "\\bar{k}")
 mat_str = mat_str.replace("k_{1} \\lambda_{10}^{2} \\left(I_{2} + I_{3}\\right)", "\\bar{k} \lambda_{10}^{2}")
 mat_tex = "\\begin{align*} \n \widetilde{G}_{%s} = \\frac{1}{2} %s \n  \end{align*} \n  \n" % (str(k),mat_str)
 all_mat_tex += mat_tex

with open("OutputGtild.txt", "w") as text_file:
 text_file.write(all_mat_tex)
