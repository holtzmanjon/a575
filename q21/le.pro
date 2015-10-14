; Jodi Berdis - Stellar Interiors Homework 2
; --------------------------------------------------------------------------------------------------------
; Write a program that solves the Lane-Emden Equation using the two boundary equations for n=3
; If y=theta and z=dtheta_dxsi, can use equations for y' and z'
; Two first order equations:
;    z' = -(y^n + (2/xsi)*z)
;    y' = z (= dtheta_dxsi)
; Initial y and z equations derived from: 
; theta_xsi = 1 - (1/6)*xsi^2 + (n/120)*xsi^4 - (n*(8*n - 5)/15120)*xsi^6
; Boundary Conditions: When xsi = 0: theta_xsi = 1, and dtheta_dxsi = z = 0
; Lane-Emden Equation: (1/(xsi^2))* d/dxsi(xsi^2 * dtheta_dxsi) = -theta^n
; --------------------------------------------------------------------------------------------------------

; Define function that returns derivatives for Lane-Emden equations
function le_diffeq,x,y,n=n
  ; given independent variable x, return N-dimensional derivatives
  ; y[0] is z, y[1] is y
  z=y[0]
  yy=y[1]
  return,[-((yy^n) + (2./x)*z),z]
end
pro rk4,x,y,h,func,_extra=extra
  ; integrate N-dimensional variable y, using stepsize h
  k1 = h * call_function(func, x, y, _extra=extra)
  k2 = h * call_function(func, x+(h/2.), y+(k1/2.), _extra=extra)
  k3 = h * call_function(func, x+(h/2.), y+(k2/2.), _extra=extra)
  k4 = h * call_function(func, x+h, y+k3, _extra=extra)
  y += (1./6.)*(k1 + (2.*k2) + (2.*k3) + k4)
end

; With functions defined, start program:
pro laneemden

; Define necessary variables:
n = 3.0                              ; polytropic index
h = 0.01                             ; step size for Runge-Kutta
xsi = 0.01                           ; set starting value of xsi to initial value of xsi
xsiarr = xsi                         ; set up blank array for xsi values

; Determine initial values of y and z (can't be exactly 0) -- comes from expansion of Lane-Emden eqn:
y=fltarr(2)
y[0] = (-(1./3.)*xsi) + (n/30.)*(xsi^3.) - (n*(8.*n - 5.)/2520.)*(xsi^5.)
y[1] = 1. - (1./6.)*(xsi^2.) + (n/120.)*(xsi^4.) - (n*(8.*n - 5.)/15120.)*(xsi^6.)
yarr = y       ; set starting value of y to initial value of y

; Do Runge-Kutta method until y < 0 (where y=0 represents the surface of the star)
WHILE y[1] GE 0 do begin
  rk4,xsi,y,h,'le_diffeq',n=n
  xsiarr = [xsiarr,xsi]
  yarr = [[yarr],[y]]
  xsi += h      ; increment xsi by the step-size value at each iteration
ENDWHILE

; ---------------------------------------------------------------------------------------------
; Part A: Calculating values for my specific polytropic index (n=3.0)

nsteps=n_elements(xarr)
xsi_1 = xsiarr[nsteps-1]         ; xsi at the surface (y=0)
z_1 = yarr[0,nsteps-1]           ; z at the surface (y=0)

G = 6.674E-8          ; gravitational constant: cm^3 grams^-1 s^-2
kb = 1.381E-16        ; boltzmann constant: erg K^-1
m_u = 1.674E-24       ; mass of a proton: grams

M = 1.989E33          ; mass of the sun: grams
R = 6.955E10          ; radius of the sun: cm
pi = 3.14159          ; value of pi
V = (4./3.)*pi*(R^3)  ; calculation of the volume of the sun

EX = 0.7              ; X: composition of hydrogen
ZE = 0.02             ; Z: composition of metals

rho_c = (-(xsi_1/3.)/(z_1)) * (M/((4./3.)*pi*(R^3.)))     ; density at the center of the sun
rho_mean = M / V                                          ; mean density of the sun
rho_c_rho_mean = rho_c / rho_mean                         ; ratio of the two density values

cap_theta = 1./(-(n+1.)*xsi_1*z_1)                        ; capital theta coefficient
mu = 4./(3.+(5.*EX)-ZE)                                   ; mean molecular weight
T_c = (cap_theta*(G*M*mu*m_u))/(kb*R)                     ; temperature at the center of the sun
Wn = 1./(4.*pi*(n+1.)*(z_1^2.))                           ; W sub-n coefficient
P_c = Wn*G*(M/R^2)^2.
Nc = (((4.*pi)^(1./n))/(n+1.))*(-(xsi_1^((n+1)/(n-1)))*z_1)^((1.-n)/n)   ; N coefficient

print, 'rho c / rho mean: ', rho_c_rho_mean
print, 'rho c: ', rho_c
print, 'cap theta: ',cap_theta
print, 'Temp c: ', T_c
print, 'Wn: ', Wn
print, 'P c: ', P_c
print, 'Wn: ', Wn
print, 'N: ', Nc

print,xsi_1, z_1

stop

;-------------------------------------------------------------------------------------------
; Part B: Plotting xsi against four 'y-axis' arrays (defined below)

yn = yarr^n                                 ; y^n or theta^n, where n is the polytropic index
yn1 = yarr^(n+1)                            ; y^n+1 or theta^n+1 --> increment value of polytrop. ind. by 1
q = ((xsiarr^2) *zarr) * (1/((xsi_1^2)*z_1))   ; q: fraction of enclosed mass at any particular radius (q=1 at surface)

!P.MULTI = [0, 1, 1]
!Y.THICK = 3
!X.THICK = 3
ps_open, '/home/users/berdis/Desktop/lane-emden'
; Plot all four lines on one single full-page plot
plot, xsiarr, yarr, xrange=[0,10], yrange=[0,1], xtitle='xsi', ytitle='theta', title='Lane-Emden Solution (Polytropic Index n = 3.0)', linestyle=0, charthick=2, charsize=1.25  ;solid   xsi vs. theta
oplot, xsiarr, yn, linestyle=1   ;dotted          xsi vs. theta^n
oplot, xsiarr, yn1, linestyle=2  ;dashed          xsi vs. theta^n+1
oplot, xsiarr, q, linestyle=5    ;long dashed     xsi vs. q (fraction of enclosed mass)

LEGEND, ['________   theta','. . . . . .   theta^n', '- - - -  theta^(n+1)', '__ __ __    q'], POSITION=[6, 0.9], charthick=2.0
ps_close

end
