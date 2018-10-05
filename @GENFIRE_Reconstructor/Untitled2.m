dim_object = dim(initialObject);
d1 = dim_object(1);
d2 = dim_object(2);
d3 = dim_object(3);

[XX,YY,ZZ] = meshgrid(np.arange(d1),np.arange(d1),np.arange(d1));
x_cen = d1/2;
y_cen = d2/2;
z_cen = d3/2;
R2 = (XX-x_cen).^2 + (YY-y_cen).^2 + (ZZ-z_cen).^2;
kernel = exp(-R2/(500^2));

print('Starting DOUGLAS RACHFORD ALGORITHM')
k = fftn(initialObject);
u = np.copy(initialObject);
for iterationNum =1:numIterations+1 %iterations are counted started from 1
    dt = 0.1 %+ 0.3*(1-np.sqrt((numIterations-iterationNum)/numIterations))
    
    if iterationNum == iterationNumsToChangeCutoff(currentCutoffNum) %update current Fourier constraint if appropriate
        constraintInd_complex = find(constraintIndicators>(obj.constraintEnforcementDelayIndicators(currentCutoffNum))&obj.measuredK~=0&obj.measuredK_mask);
        constraintInd_complex_shifted = My_iffshift3_ind(size(paddedSupport),constraintInd_complex);
        currentCutoffNum = currentCutoffNum+1;
        bestErr = 1e30;%reset best error
    end
    if enforce_positivity
        initialObject(initialObject<0) = 0; %enforce positivity
    end
    if enforce_support
        initialObject = initialObject.*paddedSupport;%enforce support
    end
    
    %take FFT of current reconstruction
    k = fftn(initialObject)
    
    %compute error
    obj.errK(iterationNum) = sum(abs(abs(k(errInd_shifted))-abs(obj.measuredK(errInd))))./sum(abs(obj.measuredK(errInd)));
    if verbose:
        print("Iteration number: {0}/{1}           error = {2:0.5f}".format(iterationNum, numIterations, errK[iterationNum-1]))
    end
    %update best object if a better one has been found
    if errK[iterationNum-1] < bestErr:
        bestErr = errK[iterationNum-1]
        outputs['reconstruction'] = initialObject
    end
    %calculate Rfree for each spatial frequency shell if necessary
    if R_freeInd_complex
        total_Rfree_error      = 0
        total_Rfree_error_norm = 0
        for shellNum in range(0, np.shape(R_freeInd_complex)[2]):
            
            tmpIndX = R_freeInd_complex[0][0][shellNum]
            tmpIndY = R_freeInd_complex[1][0][shellNum]
            tmpIndZ = R_freeInd_complex[2][0][shellNum]
            
            tmpVals = np.abs(R_freeVals_complex[shellNum])
            % Rfree_numerator                         = np.sum(abs(k[tmpIndX, tmpIndY, tmpIndZ] - tmpVals))
            Rfree_numerator                         = np.sum(abs(np.abs(k[tmpIndX, tmpIndY, tmpIndZ]) - tmpVals))
            Rfree_denominator                       = np.sum(abs(tmpVals))
            total_Rfree_error                      += Rfree_numerator
            total_Rfree_error_norm                 += Rfree_denominator
            Rfree_complex_bybin[shellNum, iterationNum-1] = Rfree_numerator / Rfree_denominator
        end
        Rfree_complex_total[iterationNum-1] = total_Rfree_error / total_Rfree_error_norm
    end
    %replace Fourier components with ones from measured data from the current set of constraints
    %k[constraintInd_complex] = measuredK[constraintInd_complex]
    k[constraintInd_complex] = dt*k[constraintInd_complex] + (1-dt)*measuredK[constraintInd_complex];
    u_K = ifftn(k);
    initialObject = 2*u_K - u;
    
    % method1
    initialObject = np.real(initialObject)
    index_neg = (initialObject<0) | (support==0);
    initialObject[index_neg]=0;
    %F_obj = np.fft.fftn(initialObject) * kernel
    %initialObject = np.real(np.fft.ifftn(F_obj))
    
    %method 2
    % object_temp = np.copy(initialObject);
    % index_re = (support) & (np.real(object_temp)>0)
    % object_temp[index_re] -= np.real(object_temp[index_re])
    % object_temp *= kernel
    % initialObject = initialObject - object_temp
    % initialObject = np.real(initialObject)
    
    u = u + initialObject - u_K
    
    %update display
    if displayFigure.DisplayFigureON:
        if iterationNum % displayFigure.displayFrequency == 0:
            if verbose:
                print("n_half_x = ", n_half_x)
                print("half_window_y = ", half_window_y)
                plot_result(initialObject,half_size,numIterations,iterationNum,errK,R_freeInd_complex,Rfree_complex_bybin,Rfree_complex_total)
                
            end %end plot_result
        end
    end
end %end methods
