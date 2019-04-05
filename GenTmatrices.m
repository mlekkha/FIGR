function GenTmatrices( numTmatrices, method, numUpstreamRegs, numDownstreamRegs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

points1 = rand(2*numTmatrices, 2);
x1 = points1(:,1);
y1 = points1(:,2);

points2 = rand(2*numTmatrices, 2);
x2 = points2(:,1);
y2 = points2(:,2);

xcoord = 0.5;
ycoord = 0.5;

Tparams = [];

for n = 1:numTmatrices
    
    if ((numDownstreamRegs==2) && (numUpstreamRegs==0))
        
        if strcmp(method, 'totRand')
            
            %From Slope-Intercept Form: y = m1x+c1
            %I define m1 = A1/B1, c1 = y1 - (A1/B1)*x1
            %After simplification, the switching boundary is: A1x - B1y + (B1*y1-A1*x1) = 0
            %Thus, T11 = A1, T12 = -B1, h1 = (B1*y1-A1*x1)
            A1 = 0;
            A1 = y1(2*n-1)-y1(2*n);
            B1 = 0;
            B1 = x1(2*n-1)-x1(2*n);
            %Define 'h1'
            h1 = 0;
            h1 = B1*y1(2*n-1)-A1*x1(2*n-1);
            
            %From Slope-Intcercept Form: y = m2x+c2
            %I define m2 = A2/B2, c2 = y2 - (A2/B2)*x2
            %After simplification, the switching boundary is: A2x - B2y + (B2*y2-A2*x2) = 0
            %Thusm T21 = A2, T22 = -B2, h2 = (B2*y2-A2*x2)
            A2 = 0;
            A2 = y2(2*n-1)-y2(2*n);
            B2 = 0;
            B2 = x2(2*n-1)-x2(2*n);
            %Define 'h2'
            h2 = 0;
            h2 = B2*y2(2*n-1)-A2*x2(2*n-1);
            
            %Define first row of Tmatrix in terms of A1 and B1
            T11 = 0;
            T11 = A1;
            T12 = 0;
            T12 = -B1;
            
            %Define first row of Tmatrix in terms of A2 and B2
            T21 = 0;
            T21 = A2;
            T22 = 0;
            T22 = -B2;
            
            Tparams = [T11,T12,h1;T21,T22,h2];
            
        elseif strcmp(method, 'goodSymmetryRand')
            
            %This method ensures that the switching boundaries intersect in
            %the middle of the phase space
            
            %From Slope-Intercept Form: y = m1x+c1
            %I define m1 = A1/B1, c1 = y1 - (A1/B1)*x1
            %After simplification, the switching boundary is: A1x - B1y + (B1*y1-A1*x1) = 0
            %Thus, T11 = A1, T12 = -B1, h1 = (B1*y1-A1*x1)
            A1 = 0;
            A1 = y1(2*n-1)-ycoord;
            B1 = 0;
            B1 = x1(2*n-1)-xcoord;
            
            %Define 'h1'
            h1 = 0;
            h1 = B1*y1(2*n-1)-A1*x1(2*n-1);
            
            %From Slope-Intcercept Form: y = m2x+c2
            %I define m2 = A2/B2, c2 = y2 - (A2/B2)*x2
            %After simplification, the switching boundary is: A2x - B2y + (B2*y2-A2*x2) = 0
            %Thusm T21 = A2, T22 = -B2, h2 = (B2*y2-A2*x2)
            A2 = 0;
            A2 = y2(2*n-1)-ycoord;
            B2 = 0;
            B2 = x2(2*n-1)-xcoord;
            %Define 'b2'
            h2 = 0;
            h2 = B2*y2(2*n-1)-A2*x2(2*n-1);
            
            %Define first row of Tmatrix in terms of m and b
            T11 = 0;
            T11 = A1;
            T12 = 0;
            T12 = -B1;
            
            %Define first row of Tmatrix in terms of m and b
            T21 = 0;
            T21 = A2;
            T22 = 0;
            T22 = -B2;
            
            Tparams = [T11,T12,h1;T21,T22,h2];
                        
        elseif strcmp(method, 'badSymmetryRand')
            
            for j = 1:numDownstreamRegs
                
                vec = [];
                
                while (size(vec,1)<1)
                    
                    testVec = rand(1,numDownstreamRegs);
                    eval = nDimCircleEval(testVec,numDownstreamRegs,0.7);
                    
                    if (eval>0)
                        
                        vec = [vec;testVec];
                        
                    else
                        
                    end
                    
                end
                
                Tparams(j,:) = [(0.5-vec(1)),(0.5-vec(2)),vec(1)*(vec(1)-0.5)+vec(2)*(vec(2)-0.5)];
                
            end
            
        elseif strcmp(method, 'autoActivationRand')
            
            while (size(Tparams,1) < 2)
                
                testMat = rand(2,2);
                
                testT = [testMat(2,2)-testMat(1,2),testMat(1,1)-testMat(2,1)];
                
                if (testT(size(Tparams,1)+1) > 0)
                    
                    Tparams = [Tparams;testT,testMat(2,1)*(testMat(1,2)-testMat(2,2))+testMat(2,2)*(testMat(2,1)-testMat(1,1))];
                    
                else
                    
                end
                
            end
            
        elseif strcmp(method, 'autoRepressionRand')
            
            while (size(Tparams,1) < 2)
                
                testMat = rand(2,2);
                
                testT = [testMat(2,2)-testMat(1,2),testMat(1,1)-testMat(2,1)];
                
                if (testT(size(Tparams,1)+1) < 0)
                    
                    Tparams = [Tparams;testT,testMat(2,1)*(testMat(1,2)-testMat(2,2))+testMat(2,2)*(testMat(2,1)-testMat(1,1))];
                    
                else
                    
                end
                
            end
            
        end
                
    elseif ((numDownstreamRegs+numUpstreamRegs)>2)
        
        if strcmp(method, 'crossProd')
            
            for j = 1:numDownstreamRegs
                
                mat = rand(numUpstreamRegs + numDownstreamRegs, numUpstreamRegs + numDownstreamRegs);
                dispMat = [];
                
                for i=1:(numDownstreamRegs+numUpstreamRegs-1)
                    
                    dispVec = mat(1,:)-mat(i+1,:);
                    dispMat = [dispMat;dispVec];
                    
                end
                
                normVec = myCrossProd(dispMat);
                Tparams(j,:) = [normVec,(-1)*dot(normVec,mat(1,:))];
                
            end
            
        elseif strcmp(method, 'badSymmetryCrossProd')
            
            for j = 1:numDownstreamRegs
                
                vec = [];
                
                while (size(vec,1)<1)
                    
                    testVec = rand(1,numUpstreamRegs + numDownstreamRegs);
                    eval = nDimCircleEval(testVec, numUpstreamRegs + numDownstreamRegs,0.7);
                    
                    if (eval>0)
                        
                        vec = [vec;testVec];
                        
                    else
                        
                    end
                    
                end
                %Here I am constructing a plane tangent to a
                %sphere of radius >= 0.7, centered in the
                %middle of the phase space.
                normVec = 2*vec; %equivalent of evaluating the gradient of the equation of a sphere at the point specified by vec
                h = dot(normVec, -vec); %this part determines the constant of the plane equation
                Tparams(j,:) = [normVec,h];
                            
            end
            
        elseif strcmp(method, 'goodSymmetryCrossProd')
            
            for j = 1:numDownstreamRegs
                
                mat = rand(numUpstreamRegs + numDownstreamRegs, numUpstreamRegs + numDownstreamRegs);
                mat(1,:) = 0.5*ones(1,numUpstreamRegs+numDownstreamRegs);
                dispMat = [];
                
                for i=1:(numDownstreamRegs+numUpstreamRegs-1)
                    
                    dispVec = mat(1,:)-mat(i+1,:);
                    dispMat = [dispMat;dispVec];
                    
                end
                
                normVec = myCrossProd(dispMat);
                Tparams(j,:) = [normVec,(-1)*dot(normVec,mat(1,:))];
                
            end
            
        elseif strcmp(method, 'autoActivationCrossProd')
            
            while (size(Tparams,1) < numDownstreamRegs)
                
                mat = [];
                mat = rand(numUpstreamRegs + numDownstreamRegs, numUpstreamRegs + numDownstreamRegs);
                dispMat = [];
                
                for i=1:(numDownstreamRegs+numUpstreamRegs-1)
                    
                    dispVec = mat(1,:)-mat(i+1,:);
                    dispMat = [dispMat;dispVec];
                    
                end
                
                normVec = [];
                normVec = myCrossProd(dispMat);
                
                if (normVec(size(Tparams,1)+1) > 0)
                    
                    Tparams = [Tparams;normVec,(-1)*dot(normVec,mat(1,:))];
                    
                else
                    
                end
                
            end
            
        elseif strcmp(method, 'autoRepressionCrossProd')
            
            while (size(Tparams,1) < numDownstreamRegs)
                
                mat = [];
                mat = rand(numUpstreamRegs + numDownstreamRegs, numUpstreamRegs + numDownstreamRegs);
                dispMat = [];
                
                for i=1:(numDownstreamRegs+numUpstreamRegs-1)
                    
                    dispVec = mat(1,:)-mat(i+1,:);
                    dispMat = [dispMat;dispVec];
                    
                end
                
                normVec = [];
                normVec = myCrossProd(dispMat);
                
                if (normVec(size(Tparams,1)+1) < 0)
                    
                    Tparams = [Tparams;normVec,(-1)*dot(normVec,mat(1,:))];
                    
                else
                    
                end
                
            end
            
        end
        
    end
     
    if numDownstreamRegs==2
        
        %Define $input as Mat
        Mat = [];
        
        if (numUpstreamRegs == 0)
            
            maternalConnects = sprintf('%f\t%f\n',[zeros(1,numDownstreamRegs)]);
            extStrengths = '\n';
            
        elseif (numUpstreamRegs == 1)
            
            maternalConnects = sprintf('%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = '\n';
            
        elseif (numUpstreamRegs == 2)
            
            maternalConnects = sprintf('%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = [sprintf('%f\n',Tparams(1,numDownstreamRegs+2)),sprintf('%f\n', Tparams(2,numDownstreamRegs+2))];
            
        elseif (numUpstreamRegs == 3)
            
            maternalConnects = sprintf('%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = [sprintf('%f\t%f\n', [Tparams(1,numDownstreamRegs+2:numDownstreamRegs+3)]),sprintf('%f\t%f\n', [Tparams(2,numDownstreamRegs+2:numDownstreamRegs+3)])];
                        
        end
        
        Mat = [sprintf(['$input\n','promoter_strengths:\n',sprintf('%f\t%f\n',log(2)/50,log(2)/50),'genetic_interconnect_matrix:\n',sprintf('%f\t%f\n',[Tparams(1,1:numDownstreamRegs)]),sprintf('%f\t%f\n',[Tparams(2,1:numDownstreamRegs)]),'external_input_strengths:\n',extStrengths,'\n','maternal_connection_strengths:\n',maternalConnects,'promoter_thresholds:\n',sprintf('%f\t%f\n',[Tparams(:,end)]),'diffusion_parameter(s):\n',sprintf('%f\n',0),'protein_half_lives:\n',sprintf('%f\t%f\n',50,50),'translational_transcriptional_delays:\n',sprintf('%f\t%f\n',6,6),'$$'])];
        
        filename = ['~/genecircuits/Tmatrices/tmat' num2str(n)];
        dlmwrite(filename, Mat, 'delimiter',  '');
        %cd genecircuits/Tmatrices
        %command = sprintf('./cat_tmat_seedn16-input.bash %d', n);
        %system(command);
        %cd ..
        %cd ..
        
        
        
        
        
        
        
        
        
        
        
        
    elseif numDownstreamRegs==3
        
        Mat = [];
        
        if (numUpstreamRegs == 0)
            
            maternalConnects = sprintf('%f\t%f\t%f\n',[zeros(1,numDownstreamRegs)]);
            extStrengths = '\n';
            
        elseif (numUpstreamRegs == 1)
            
            maternalConnects = sprintf('%f\t%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = '\n';
            
        elseif (numUpstreamRegs == 2)
            
            maternalConnects = sprintf('%f\t%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = [sprintf('%f\n', [Tparams(1,numDownstreamRegs+2)]), sprintf('%f\n', [Tparams(2,numDownstreamRegs+2)]), sprintf('%f\n', [Tparams(3,numDownstreamRegs+2)])];
            
        elseif (numUpstreamRegs == 3)
            
            maternalConnects = sprintf('%f\t%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = [sprintf('%f\t%f\n', [Tparams(1,numDownstreamRegs+2:numDownstreamRegs+3)]), sprintf('%f\t%f\n', [Tparams(2,numDownstreamRegs+2:numDownstreamRegs+3)]), sprintf('%f\t%f\n', [Tparams(3,numDownstreamRegs+2:numDownstreamRegs+3)])];
            
        end
        
        Mat = [sprintf(['$input\n','promoter_strengths:\n',sprintf('%f\t%f\t%f\n',log(2)/50,log(2)/50,log(2)/50),'genetic_interconnect_matrix:\n',sprintf('%f\t%f\t%f\n',Tparams(1,1:numDownstreamRegs)),sprintf('%f\t%f\t%f\n',Tparams(2,1:numDownstreamRegs)),sprintf('%f\t%f\t%f\n',Tparams(3,1:numDownstreamRegs)),'external_input_strengths:\n',extStrengths,'\n','maternal_connection_strengths:\n',maternalConnects,'promoter_thresholds:\n',sprintf('%f\t%f\t%f\n',[Tparams(:,end)]),'diffusion_parameter(s):\n',sprintf('%f\n',0),'protein_half_lives:\n',sprintf('%f\t%f\t%f\n',50,50,50),'translational_transcriptional_delays:\n',sprintf('%f\t%f\t%f\n',6,6,6),'$$'])];
        filename = ['~/genecircuits/Tmatrices/tmat' num2str(n)];
        dlmwrite(filename, Mat, 'delimiter',  '');
        cd genecircuits/Tmatrices
        command = sprintf('./cat_tmat_seedn16-input.bash %d', n);
        system(command);
        cd ..
        cd ..
        
    elseif numDownstreamRegs==4
        
        Mat = [];
        
        if (numUpstreamRegs == 0)
            
            maternalConnects = sprintf('%f\t%f\t%f\t%f\n',[zeros(1,numDownstreamRegs)]);
            extStrengths = '\n';
            
        elseif (numUpstreamRegs == 1)
            
            maternalConnects = sprintf('%f\t%f\t%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = '\n';
            
        elseif (numUpstreamRegs == 2)
            
            maternalConnects = sprintf('%f\t%f\t%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = [sprintf('%f\n', [Tparams(1,numDownstreamRegs+2)]), sprintf('%f\n', [Tparams(2,numDownstreamRegs+2)]), sprintf('%f\n', [Tparams(3,numDownstreamRegs+2)]),sprintf('%f\n', [Tparams(4,numDownstreamRegs+2)])];
            
        elseif (numUpstreamRegs == 3)
            
            maternalConnects = sprintf('%f\t%f\t%f\n',[Tparams(:,numDownstreamRegs+1)]);
            extStrengths = [sprintf('%f\t%f\n', [Tparams(1,numDownstreamRegs+2:numDownstreamRegs+3)]), sprintf('%f\t%f\n', [Tparams(2,numDownstreamRegs+2:numDownstreamRegs+3)]), sprintf('%f\t%f\n', [Tparams(3,numDownstreamRegs+2:numDownstreamRegs+3)]),sprintf('%f\t%f\n', [Tparams(4,numDownstreamRegs+2:numDownstreamRegs+3)])];
            
        end
        
        Mat = [sprintf(['$input\n','promoter_strengths:\n',sprintf('%f\t%f\t%f\t%f\n',log(2)/50,log(2)/50,log(2)/50,log(2)/50),'genetic_interconnect_matrix:\n',sprintf('%f\t%f\t%f\t%f\n',Tratios(1,1), Tratios(1,2),Tratios(1,3),Tratios(1,4)),sprintf('%f\t%f\t%f\t%f\n',Tratios(2,1),Tratios(2,2),Tratios(2,3),Tratios(2,4)),sprintf('%f\t%f\t%f\t%f\n',Tratios(3,1),Tratios(3,2),Tratios(3,3),Tratios(3,4)),sprintf('%f\t%f\t%f\t%f\n',Tratios(4,1), Tratios(4,2),Tratios(4,3),Tratios(4,4)),'external_input_strengths:\n',extStrengths,'\n','maternal_connection_strengths:\n',maternalConnects,'promoter_thresholds:\n',sprintf('%f\t%f\t%f\t%f\n',[Tparams(:,end)]),'diffusion_parameter(s):\n',sprintf('%f\n',0),'protein_half_lives:\n',sprintf('%f\t%f\t%f\t%f\n',50,50,50,50),'translational_transcriptional_delays:\n',sprintf('%f\t%f\t%f\t%f\n',6,6,6,6),'$$'])];
        
        filename = ['~/genecircuits/Tmatrices/tmat' num2str(n)];
        dlmwrite(filename, Mat, 'delimiter',  '');
        cd genecircuits/Tmatrices
        command = sprintf('./cat_tmat_seedn16-input.bash %d', n);
        system(command);
        cd ..
        cd ..
        
    end
    fprintf ('asdfasdf');

end

end

