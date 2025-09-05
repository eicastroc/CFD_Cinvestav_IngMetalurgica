using LinearAlgebra

"""
heatTransient()

Método para resolver para n pasos de tiempo, Δt, dados, 
el sistema de ecuaciones que describe la aplicación en régimen transitorio: 
	
		[A][x] = [b], 
donde
	[A] es la matriz de coeficientes, 
	[x] es un vector que contiene las temperaturas en los volúmenes de control,
	[b] es el vector constante al lado derecho de la ecuacion.

Parámetros:
	- X (float): Tamaño del dominio, en [m].
	- nCells (int): Número de celdas.
	- Δt (float): tamaño de paso de tiempo, en [s].
	- n_Δt (int): número de pasos de tiempo a efectuar en el cálculo.
	- Told (Vector, float): perfil de temperatura inicial, de tamaño nCells, en [K].
	- λ (float): conductividad térmica, en [W.m-1.K-1].
	- ρ (float): densidad, en [kg.m-3].
	- c (float): capacidad calorífica, en [J.kg-1.K-1].
	- htc (float): coeficiente de transferencia de calor por convección, en [W.m-2.K-1]
	- TfarW (float): temperatura "lejos" en la frontera izquierda, en [K].
	- TfarE (float): temperatura "lejos" en la frontera derecha, en [K].
	- niters (int): número máximo de iteraciones.
	- tol (float): tolerancia de convergencia del método numérico.

Regresa:
	- Tnew (vector, float): perfil de temperatura en régimen transitorio, de tamaño nCells, en [K].
	- error (vector, float): error de convergencia, de tamaño del número de pasos de tiempo.
	- bilan (vector, float): balance de energía al lograr convergencia en cada paso de tiempo.

"""
function heatTransient(X, nCells, Δt, n_Δt, Told, λ, ρ, c, htc, TfarW, TfarE, niters, tol=1e-15)

	error = ones(Float64, n_Δt, 2) # para cada Δt: numero de its., convergencia.
	bilan = ones(Float64, n_Δt) # para cada Δt: error en balance de energia

	# Inicializacion de variables
	Tstar = copy(Told) # temperatura iterada hasta convergencia
	Tnew = copy(Told)  # nueva temperatura calculada
	iteracion = 1
	err = 1
	
	for n in 1:n_Δt
		# Valores iniciales para el paso de tiempo
		Tstar = copy(Told)
		Tnew = copy(Told) 
		
		###### Método de avance en un paso de tiempo, Δt
		for i=1:niters
			# Construccion del sistema de ecuaciones
			A = constructLHS_trans(X, nCells, Tstar, Δt, λ, ρ, c, htc)
			b = constructRHS_trans(X, nCells, Told, Tstar, Δt, λ, ρ, c, htc, TfarW, TfarE)
			# Resolucion del sistema
			Tnew = A\b
	
			# Estimacion del error numerico
			err = abs(maximum((Tnew .- Tstar) ./ Told))

			# Actualizacion de la temperatura iterada
			Tstar = Tnew
	
			if err <= tol
				iteracion = i
				if err < 1e-16
					err = 1e-16 # para poder graficar
				end
				break
			end
		end
		
		# Registro del error numérico
		error[n,1] = iteracion
		error[n,2] = err
		
		###### Balance de energía
	 	# Calculo temperaturas en la frontera
		T_bw = tempBound(X, nCells, λ, htc, TfarW, Tnew[begin])
		T_be = tempBound(X, nCells, λ, htc, TfarE, Tnew[end])
				
		# Estimacion de error balance de energia
		ϕ_e, ϕ_w, Gen, Acc = transientHeatBalance(X, nCells, Δt, λ, ρ, c, Told, Tnew, T_bw, T_be)
		bil = abs(ϕ_e + ϕ_w + Gen - Acc)
		if bil < 1e-16
			bil = 1e-16 # para poder graficar
		end
		bilan[n] = bil
	
		# Actualizacion de temperatura para próximo paso de tiempo, Δt
		Told = copy(Tnew)

	end
	
	return Tnew, error, bilan

end


"""
Sc_()

Función para calcular la parte constante, Sc, 
de la ecuación de un término fuente linealizado, 
de la forma: S = Sc + Sp * T

Parámetros:
	- T (float): valor de T* 

Regresa:
	- Sc (float): valor de Sc

"""
function Sc_(T)
	return 0.0
end


"""
Sp_()

Función para calcular la pendiente, Sp, 
de la ecuación de un término fuente linealizado, 
de la forma: S = Sc + Sp * T

Parámetros:
	- T (float): valor de T* 

Regresa:
	- Sp (float): valor de Sp

"""
function Sp_(T)
	return 0.0
end


"""
constructLHS_trans()

Construcción de la matriz de coeficientes [A], al lado izquierdo 
del sistema de ecuaciones para régimen transitorio.

Parámetros
	- X (float): Tamaño del dominio, en [m].
	- nCells (int): Número de celdas.
	- Tstar (Vector, float): perfil de temperaturas iteradas, de tamaño nCells, en [K].
	- Δt (float): paso de tiempo, en [s].
	- λ (float): conductividad térmica, en [W.m-1.K-1].
	- ρ (float): densidad, en [kg.m-3].
	- c (float): capacidad calorífica, en [J.kg-1.K-1].
	- htc (float): coeficiente de transferencia de calor por convección, en [W.m-2.K-1]

Regresa:
	- Tridiagonal(float): matriz tridiagonal con coeficientes del sistema de ecuaciones.
"""
function constructLHS_trans(X, nCells, Tstar, Δt, λ, ρ, c, htc)

	# Discretizacion espacial
	Δx = X/nCells 	# [m]
	A = 1 			# [m2]
	ΔV = Δx * A 	# [m3]


	# Coeficientes ap0:
	ap0 = ρ * c * ΔV / Δt
	
	# Coeficientes aE, aW:
	aE = λ * A / Δx
	aW= λ * A / Δx

	# Resistencia a las fronteras por conveccion + difusion
	Rth = (htc * (λ/(Δx/2))) / (htc + (λ/(Δx/2))) 
		
	# diagonal inferior de la matriz A (coeficientes aW)
	dl = @. -aW * ones(Float64, nCells-1)
	# diagonal superior de la matriz A (coeficientes aE)
	du = @. -aE * ones(Float64, nCells-1)
	# diagonal de la matriz A (coeficientes aP)
	d = zeros(Float64, nCells)
	d[1] = ap0 + aE + Rth * A - Sp_(Tstar[1]) * ΔV
	for i=2:nCells-1
		d[i] = ap0 + aE + aW - Sp_(Tstar[i]) * ΔV
	end
	d[nCells] = ap0 + aE + Rth * A - Sp_(Tstar[nCells]) * ΔV

	return Tridiagonal(dl, d, du)
end


"""
constructRHS_trans()

Construcción del vector constante [b], al lado derecho 
del sistema de ecuaciones en régimen estacionario.

Parámetros
	- X (float): Tamaño del dominio, en [m].
	- nCells (int): Número de celdas.
	- Told (Vector, float): perfil de temperatura inicial, de tamaño nCells, en [K].
	- Tstar (Vector, float): perfil de temperaturas iteradas, de tamaño nCells, en [K].
	- Δt: paso de tiempo, en [s].
	- λ (float): conductividad térmica, en [W.m-1.K-1].
	- ρ (float): densidad, en [kg.m-3].
	- c (float): capacidad calorífica, en [J.kg-1.K-1].
	- htc (float): coeficiente de transferencia de calor por convección, en [W.m-2.K-1]
	- Tw (float): temperatura "lejos" en la frontera izquierda, en [K].
	- Te (float): temperatura "lejos" en la frontera derecha, en [K].

Regresa:
	- Vector(float): vector con coeficientes del sistema de ecuaciones.
"""
function constructRHS_trans(X, nCells, Told, Tstar, Δt, λ, ρ, c, htc, Tw, Te)

	# Discretizacion espacial
	Δx = X/nCells 	# [m]
	A = 1 			# [m2]
	ΔV = Δx * A 	# [m3]

	# Coeficientes ap0:
	ap0 = ρ * c * ΔV / Δt
	
	# Resistencia a las fronteras por conveccion + difusion
	Rth = (htc * (λ/(Δx/2))) / (htc + (λ/(Δx/2))) 
	
	RHS = zeros(Float64, nCells)

	RHS[1] = ap0 * Told[1] + Rth*A*Tw + Sc_(Tstar[1]) * ΔV
	for i=2:nCells-1
		RHS[i] = ap0 * Told[i] + Sc_(Tstar[i]) * ΔV
	end
	RHS[nCells] = ap0 * Told[nCells] + Rth*A*Te + Sc_(Tstar[nCells]) * ΔV

	return Vector(RHS)
end


"""
tempBound()

Calculo de la temperatura en la frontera para una 
condición a la frontera tipo Cauchy, 
mediante un método de resistencias térmicas equivalentes.

Parámetros
	- X (float): Tamaño del dominio, en [m].
	- nCells (int): Número de celdas.
	- λ (float): conductividad térmica, en [W.m-1.K-1].
	- htc (float): coeficiente de transferencia de calor por convección, en [W.m-2.K-1]
	- Tinf (float): temperatura "lejos" en la frontera, en [K].
	- Tcell (float): temperatura en la celda en la frontera, en [K].

Regresa:
	- (float): temperatura en la frontera, en [K]
"""
function tempBound(X, nCells, λ, htc, Tinf, Tcell)

	# Discretizacion espacial
	Δx = X/nCells 	# [m]

	# Resistencias termicas
	rConv = htc  		# W/(m2 K)
	rDiff = λ / (Δx/2) 	# W/(m2 K)
	
	return (rConv*Tinf + rDiff*Tcell) / (rConv + rDiff)
end


"""
transientHeatBalance()

Balance de energía para el régimen transitorio, en [W]:

	Acumulación = Entradas - Salidas + Generación

Parámetros
	- X (float): Tamaño del dominio, en [m].
	- nCells (int): Número de celdas.
	- Δt: paso de tiempo, en [s].
	- λ (float): conductividad térmica, en [W.m-1.K-1].
	- ρ (float): densidad, en [kg.m-3].
	- c (float): capacidad calorífica, en [J.kg-1.K-1].
	- Told (vector, float): perfil de temperatura paso de tiempo previo, de tamaño nCells, en [K].
	- Tnew (vector, float): perfil de temperatura, de tamaño nCells, en [K].
	- Twest (float): temperatura en la frontera izquierda, en [K].
	- Teast (float): temperatura en la frontera derecha, en [K].

Regresa:
	- (float): temperatura en la frontera, en [K]
"""
function transientHeatBalance(X, nCells, Δt, λ, ρ, c, Told, Tnew, Twest, Teast)

	# Discretizacion espacial
	Δx = X/nCells 	# [m]
	A = 1 			# [m2]
	ΔV = Δx * A 	# [m3]

	# Flujo en la frontera oeste
	ϕ_w = -λ/(Δx/2) * A * (Tnew[begin] - Twest)
	# Flujo en la frontera este
	ϕ_e = λ/(Δx/2) * A * (Teast - Tnew[end])
	# Generacion
	Gen = sum(@. ΔV * (Sc_(Tnew) + Sp_(Tnew) * Tnew))
	# Acumulacion
	Acc = sum(@. ρ * c * ΔV / Δt * (Tnew - Told))
	
	return ϕ_w, ϕ_e, Gen, Acc
end


"""
CtoK()

Conversión de unidades de temperatura de Celcius a Kelvin.

Parámetros
	- T (float): Temperatura, en [C].

Regresa:
	- (float): Temperatura, en [K].

"""
function CtoK(T)
	# convert temperature: Celcius to Kelvin
	return @. T + 273.15
end


"""
KtoC()

Conversión de unidades de temperatura de Kelvin a Celcius

Parámetros
	- T (float): Temperatura, en [K].

Regresa:
	- (float): Temperatura, en [C].
"""
function KtoC(T)
	# convert temperature: Kelvin to Celcius
	return @. T - 273.15
end


println("heatTransient.jl has been imported")