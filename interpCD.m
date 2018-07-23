function [CDi] = interpCD(Mi,Ai,CLi)

%1-h~m,2-T~deg C,3-g~m/s2,4-P~N/m2,5-rho~kg/m3,6-mu~Ns/m2
%http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
% DATA=[-1000	 21.50	9.810	11.39       13.47   	1.821
%         0	 15.00	9.807	10.13       12.25   	1.789
%        1000	  8.50	9.804	8.988       11.12   	1.758
%        2000	  2.00	9.801	7.950       10.07   	1.726
%        3000	 -4.49	9.797	7.012       9.093   	1.694
%        4000	-10.98	9.794	6.166       8.194   	1.661
%        5000	-17.47	9.791	5.405       7.364   	1.628
%        6000	-23.96	9.788	4.722       6.601   	1.595
%        7000	-30.45	9.785	4.111       5.900   	1.561
%        8000	-36.94	9.782	3.565       5.258   	1.527
%        9000	-43.42	9.779	3.080       4.671   	1.493
%       10000	-49.90	9.776	2.650       4.135   	1.458
%       15000	-56.50	9.761	1.211       1.948       1.422
%       20000	-56.50	9.745	0.5529  	0.8891      1.422
%       25000	-51.60	9.730	0.2549  	0.4008      1.448
%       30000	-46.64	9.715	0.1197  	0.1841      1.475
%       40000	-22.80	9.684	0.0287      0.03996 	1.601
%       50000	-25.00	9.654	0.007978	0.01027     1.704
%       60000	-26.13	9.624	0.002196	0.003097	1.584
%       70000	-53.57	9.594	0.00052 	0.0008283	1.438
%       80000	-74.51	9.564	0.00011 	0.0001846	1.321];
    
%M Matrix
Mp(:,:,1) =[0.3477                   0.54385                    0.74                    0.78                      0.82
            0.3477                   0.54385                    0.74                    0.78                      0.82
            0.3477                   0.54385                    0.74                    0.78                      0.82];


Mp(:,:,2) =[0.3477                   0.54385                    0.74                    0.78                      0.82
            0.3477                   0.54385                    0.74                    0.78                      0.82
            0.3477                   0.54385                    0.74                    0.78                      0.82];

Mp(:,:,3) =[0.3477                   0.54385                    0.74                    0.78                      0.82
            0.3477                   0.54385                    0.74                    0.78                      0.82
            0.3477                   0.54385                    0.74                    0.78                      0.82];
           
%A Matrix
Ap(:,:,1) =[6114.29534726905          16034.2211732974          25954.1469993257        27977.0734996629          30000
            6114.29534726905          16034.2211732974          25954.1469993257        27977.0734996629          30000
            6114.29534726905          16034.2211732974          25954.1469993257        27977.0734996629          30000];

Ap(:,:,2) =[9731.74006493887          21354.4067823009          32977.0734996628        33988.5367498314          35000
            9731.74006493887          21354.4067823009          32977.0734996628        33988.5367498314          35000
            9731.74006493887          21354.4067823009          32977.0734996628        33988.5367498314          35000];

Ap(:,:,3) =[13349.1847826087          26674.5923913043                     40000        40000                     40000
            13349.1847826087          26674.5923913043                     40000        40000                     40000
            13349.1847826087          26674.5923913043                     40000        40000                     40000];
     

% M=[0.3477                   0.54385                    0.74                    0.78                      0.82];
% A=[30000                    35000                      40000   ];
% CL(:,:) =[1.29168242029875         0.778056284613185         0.638479022607684        0.628480272690572         0.622835974430675
%             1.73884387443703         1.04740792566795          0.859511068642657        0.846050898579046         0.838452627285359
%             2.18600532857531         1.31675956672271          1.08054311467763         1.06362152446752          1.05406928014004];
% figure(1)
% % surf(M,A,CL)
% V = interp2(M,A,CL,Mi,Ai,'spline') 
%CL Matrix                     
CL(:,:,1) =[1.29168242029875         0.778056284613185         0.638479022607684        0.628480272690572         0.622835974430675
            1.73884387443703         1.04740792566795          0.859511068642657        0.846050898579046         0.838452627285359
            2.18600532857531         1.31675956672271          1.08054311467763         1.06362152446752          1.05406928014004];
        
CL(:,:,2) =[1.48304775843925         0.969741160806569         0.876949323180298        0.827555427779135         0.785414682112987
            1.9964570777879          1.30545128644544          1.18053627953772         1.11404294409925          1.05731369214517
            2.50986639713656         1.6411614120843           1.48412323589514         1.40053046041937          1.32921270217734];
 
CL(:,:,3) =[1.70905177605013         1.2201696876878           1.22498201893316         1.10256435497666          0.997620692397083
            2.30070036193042         1.64257448569979          1.64905277523617         1.48425591669844          1.34298230178365
            2.89234894781072         2.06497928371178          2.07312353153919         1.86594747842022          1.68834391117023];

%CD Matrix
CD(:,:,1) =[1.29168242029875         0.778056284613185         0.638479022607684        0.628480272690572         0.622835974430675
            1.73884387443703         1.04740792566795          0.859511068642657        0.846050898579046         0.838452627285359
            2.18600532857531         1.31675956672271          1.08054311467763         1.06362152446752          1.05406928014004];
        
CD(:,:,2) =[1.48304775843925         0.969741160806569         0.876949323180298        0.827555427779135         0.785414682112987
            1.9964570777879          1.30545128644544          1.18053627953772         1.11404294409925          1.05731369214517
            2.50986639713656         1.6411614120843           1.48412323589514         1.40053046041937          1.32921270217734];
 
CD(:,:,3) =[1.70905177605013         1.2201696876878           1.22498201893316         1.10256435497666          0.997620692397083
            2.30070036193042         1.64257448569979          1.64905277523617         1.48425591669844          1.34298230178365
            2.89234894781072         2.06497928371178          2.07312353153919         1.86594747842022          1.68834391117023];

Mp1=reshape(Mp,numel(Mp),1);   
Ap1=reshape(Ap,numel(Ap),1);   
CL1=reshape(CL,numel(CL),1);   
CD1=reshape(CD,numel(CD),1);   

% figure
% hold on
% for i=1:size(Mp,1)
%   plot3(reshape(Mp(i,:,:),size(Mp,2)*size(Mp,3),1),reshape(Ap(i,:,:),size(Mp,2)*size(Mp,3),1),reshape(CL(i,:,:),size(Mp,2)*size(Mp,3),1),'r-');
% end
% for i=1:size(Mp,2)
%   plot3(reshape(Mp(:,i,:),size(Mp,1)*size(Mp,3),1),reshape(Ap(:,i,:),size(Mp,1)*size(Mp,3),1),reshape(CL(:,i,:),size(Mp,1)*size(Mp,3),1),'b-');
% end
% for i=1:size(Mp,3)
%   plot3(reshape(Mp(:,:,i),size(Mp,1)*size(Mp,2),1),reshape(Ap(:,:,i),size(Mp,1)*size(Mp,2),1),reshape(CL(:,:,i),size(Mp,1)*size(Mp,2),1),'k-');
% end

% Vint = scatteredInterpolant(Mp1,Ap1,CL1,CD1);
% V = Vint(Mi,Ai,CLi) 
plot3(Mp1,Ap1,CL1,'k*');

% V = griddata(Mp1,Ap1,CL1,CD1,Mi,Ai,CLi,'linear') 
end