module XRFast

export nnls, subX, lasso, lasso_sparse, approximate_ksvd, dictionary_learning, dictionary_learning_1, LLS, INV_LLS, Savitsky_Golay, SNIP, initial_dictionary, approx_nn_dictionary_learning, full_nn_dictionary_learning
 
using Lasso
using LinearAlgebra
using DSP
using ProgressMeter
using SparseArrays
using NonNegLeastSquares

#subX function randomly selects 5% of spectra from original data cube and saves as a new array
#Spectral dimensions in y (dim 1), pixels in x (dim 2)
function subX(a)
    channels = size(a,2)
    pixels = size(a,1)
    percent = trunc(Int64, (size(a,1))*0.05) #select 5% of pixels to form 
    subX =zeros(0) #empty array
    for i in 1:percent #loop to create array of randomly selected pixels. Stored in an array called subX
        n = size(a,1)
        idx = rand(1:n)
        append!(subX,a[idx,:])
    end
    return reshape(subX,channels,percent)
end

function initial_dictionary(a;
    Atoms::Int = default_Atoms)
 
    default_Atoms = 40
 
    channels = size(a,2)
    pixels = size(a,1)
    Prc = Atoms/pixels

    percent = trunc(Int64, (size(a,1))*Prc)  
    subX =zeros(0) #empty array
    for i in 1:percent 
        n = size(a,1)
        idx = rand(1:n)
        append!(subX,a[idx,:])
    end
    return reshape(subX,channels,percent)
end

function lasso_sparse(subx, D, Threshold)
    #p = size(subx,2) # numb of pixels in image
    w = zeros(0)
    for i in 1:size(subx,2)
        q = subx[:,i]
        W = fit(LassoPath, D, q;intercept=false, cd_tol=1e-5)
        J = Matrix(coef(W))
        append!(w,J[:,size(J,2)])
    end
    W = reshape(w,size(D,2),size(subx,2))
    T = [broadcast(abs, x) > Threshold ? x : 0 for x in W] #Threshold parameter to make sparse, necessary?
    T =[val*1.0 for val in T]
    return T
end

function lasso(subx, D, cd_tolerance)

    #p = size(subx,2) # numb of pixels in image
    w = zeros(0)
    for i in 1:size(subx,2)
        q = subx[:,i]
        W = fit(LassoPath, D, q;intercept=false, cd_tol=cd_tolerance)
        J = Matrix(coef(W))
        append!(w,J[:,size(J,2)])
    end
    W = reshape(w,size(D,2),size(subx,2))
    return W
end


function approximate_ksvd(Y::AbstractMatrix, D::AbstractMatrix, W::AbstractMatrix)
    R = Y - D*W
    for k in 1:size(D,2)
        I = findall( x->(x > 0), W[k,:])
        Ri = R[:,I] + D[:,k]*W[k,I']
        dk = Ri * W[k,I]
        dk = dk/sqrt(dk'*dk)  # normalize
        D[:,k] = dk
        W[k,I] = dk'*Ri;
        R[:,I] = Ri - D[:,k]*W[k,I']
    end
    return D
end



function nn_approximate_ksvd(Y::AbstractMatrix, D::AbstractMatrix, W::AbstractMatrix)
    R = Y - D*W
    for k in 1:size(D,2)
        I = findall( x->(x > 0), W[k,:])
	    #R[:,I] = [x > 0 ? x : 0 for x in R[:,I]]
        Ri = R[:,I] + D[:,k]*W[k,I']
        dk = Ri * W[k,I]
        dk = dk/sqrt(dk'*dk)  # normalize
	    @. dk[isnan(dk)] = 0
        D[:,k] = dk
        D[:,k] = [x > 0 ? x : 0 for x in D[:,k]]
	    W[k,I] = D[:,k]'*Ri;
        R[:,I] = Ri - D[:,k]*W[k,I']
    end
    return D
end

function nn_full_ksvd(Y::AbstractMatrix, D::AbstractMatrix, W::AbstractMatrix)
    R = Y - D*W
    for k in 1:size(D,2)
        I = findall( x->(x > 0), W[k,:])
        Ri = R[:,I] + D[:,k]*W[k,I']
        U,S,V = svd(Ri)
        #@. U[isnan(U)] = 0
        D[:,k] = U[:,1]
        D[:,k] = [x > 0 ? x : 0 for x in D[:,k]]
        W[k,I] = V'*S
        R[:,I] = Ri - D[:,k]*W[k,I']
    end
    return D
end

function nnls(subx, D)
    w = zeros(0)
    for i in 1:size(subx,2)
        q = subx[:,i]
        J = nonneg_lsq(D,subx[:,i];alg=:nnls)
        append!(w,J[:,size(J,2)])
    end
    W = reshape(w,size(D,2),size(subx,2))
    return W
end


function dictionary_learning(D::AbstractMatrix, L::AbstractMatrix; 
        noIt_Sub::Int = default_It_Sub, 
        noIt_KSVD::Int = default_It_KSVD,
        Threshold::Float64 = default_Thresh)
    
    default_It_Sub = 10
    default_It_KSVD = 1
    default_Thresh = 0.25
    
    p = Progress(noIt_Sub, 1, "Optimizing Dictionary...")
    for i in 1:noIt_Sub
        y = subX(L)
        for j in 1:noIt_KSVD
            W = lasso_sparse(y,D,Threshold)
            D = approximate_ksvd(y, D, W)
        end
    next!(p)
    end
    return(D)
end

function dictionary_learning_1(D::AbstractMatrix, L::AbstractMatrix; 
        noIt_Sub::Int = default_It_Sub, 
        noIt_KSVD::Int = default_It_KSVD,
        cd_tolerance = default_It_KSVD)
    
    default_It_Sub = 10
    default_It_KSVD = 1
    default_cd_tolerance = 1e-6
    
    p = Progress(noIt_Sub, 1, "Optimizing Dictionary...")
    for i in 1:noIt_Sub
        y = subX(L)
        for j in 1:noIt_KSVD
            W = lasso(y,D,cd_tolerance)
            D = approximate_ksvd(y, D, W)
        end
    next!(p)
    end
    return(D)
end


function approx_nn_dictionary_learning(D::AbstractMatrix, L::AbstractMatrix; 
        noIt_Sub::Int = default_It_Sub, 
        noIt_KSVD::Int = default_It_KSVD)
    
    default_It_Sub = 10
    default_It_KSVD = 1
    
    p = Progress(noIt_Sub, 1, "Optimizing Dictionary...")
    for i in 1:noIt_Sub
        y = subX(L)
        for j in 1:noIt_KSVD
            W = nnls(y,D)
            D = nn_approximate_ksvd(y, D, W)
        end
    next!(p)
    end
    return(D)
end


function full_nn_dictionary_learning(D::AbstractMatrix, L::AbstractMatrix; 
        noIt_Sub::Int = default_It_Sub, 
        noIt_KSVD::Int = default_It_KSVD)
    
    default_It_Sub = 10
    default_It_KSVD = 1
    
    p = Progress(noIt_Sub, 1, "Optimizing Dictionary...")
    for i in 1:noIt_Sub
        y = subX(L)
        for j in 1:noIt_KSVD
            W = nnls(y,D)
            D = nn_full_ksvd(y, D, W)
        end
    next!(p)
    end
    return(D)
end

function LLS(a)
    channels = size(a,2)
    pixels = size(a,1)
    LLS =zeros(0) #empty array
    for i in 1:size(a,1)
        T = [x > 0 ? x : 0 for x in a[i,:]] 
        dr = [log(log(sqrt(val+1)+1)+1) for val in T] 
        append!(LLS,dr)
    end
    return reshape(LLS,channels,pixels)
end

function INV_LLS(a)
    channels = size(a,2)
    pixels = size(a,1)
    INV_LLS =zeros(0) #empty array
    for i in 1:size(a,1)
        df = [(exp(exp(val) - 1) - 1)^2 - 1 for val in a[i,:]]
        append!(INV_LLS,df)
    end
    return reshape(INV_LLS,channels,pixels)
end

function Savitsky_Golay(d::AbstractMatrix;
        windowSize::Int = default_windowSize, 
        polyOrder::Int = default_polyOrder, 
        deriv::Int = default_deriv)
    pixels = size(d,1) #number of pixels in image
    default_windowSize = 33 #must be odd
    default_polyOrder = 2 #must be less than the window size
    default_deriv = 0
    for i in 1:pixels
        x = d[i,:] 
        isodd(windowSize) || throw("Window size must be an odd integer.")
        polyOrder < windowSize || throw("Polynomial order must me less than window size.")
        halfWindow = Int( ceil((windowSize-1)/2) )
        S = zeros.(windowSize, polyOrder+1)
        for ct = 0:polyOrder
        S[:,ct+1] = (-halfWindow:halfWindow).^(ct)
        end
        G = S * pinv(S' * S)
        filterCoeffs = G[:, deriv+1] * factorial(deriv)
        paddedX = [x[1]*ones(halfWindow); x; x[end]*ones(halfWindow)]
        y = conv(filterCoeffs[end:-1:1], paddedX)
        xj = y[2*halfWindow+1:end-2*halfWindow]
        d[i,:] = [yal > 0 ? yal : 0 for yal in xj]
    end   
    return(d)
end

function SNIP(d::AbstractMatrix)
    M = 36 # distance of M
    iter = 15; # number of iterations  

    Mb = 14 # distance of M
    iterb = 5

    Mc = 7 # distance of M
    iterc = 4

    Md = 3 # distance of M
    iterd = 2
    
    for j in 1:iter, i in M+1:4096-M
            M1 = d[:,i-M]
            M2 = d[:,i+M]
            J = [(x+y)/2 for (x,y) in zip(M1, M2)]
            d[:,i] = [x < y ? x : y for (x, y) in zip(J, d[:,i])]
    end
    for j in 1:iterb, i in Mb+1:4096-Mb
            M1 = d[:,i-Mb]
            M2 = d[:,i+Mb]
            J = [(x+y)/2 for (x,y) in zip(M1, M2)]
            d[:,i] = [x < y ? x : y for (x, y) in zip(J, d[:,i])]
    end
    for j in 1:iterc, i in Mc+1:4096-Mc
            M1 = d[:,i-Mc]
            M2 = d[:,i+Mc]
            J = [(x+y)/2 for (x,y) in zip(M1, M2)]
            d[:,i] = [x < y ? x : y for (x, y) in zip(J, d[:,i])]
    end
    for j in 1:iterd, i in Md+1:4096-Md
            M1 = d[:,i-Md]
            M2 = d[:,i+Md]
            J = [(x+y)/2 for (x,y) in zip(M1, M2)]
            d[:,i] = [x < y ? x : y for (x, y) in zip(J, d[:,i])]
    end
    return(d)
end



end # module
