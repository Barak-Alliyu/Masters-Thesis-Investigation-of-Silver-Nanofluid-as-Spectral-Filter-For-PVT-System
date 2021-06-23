function [Tr,Rel,Abs]=select_filter_M(glass)

    if strcmp(glass,'anti_reflective_thin_glass_3.8mm')==1
      Tr=xlsread('data_M', 'optical properties', 'H3:H2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'I3:I2004')/100; 
      Abs=xlsread('data_M', 'optical properties', 'J3:J2004')/100; 
    end

    if strcmp(glass,'anti_reflective_thin_glass_0.7mm')==1
      Tr=xlsread('data_M', 'optical properties', 'M3:M2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'N3:N2004')/100; 
      Abs=xlsread('data_M', 'optical properties', 'O3:O2004')/100; 
    end

    if strcmp(glass,'pure_water_10mm')==1
      Tr=xlsread('data_M', 'optical properties', 'C3:C2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'D3:D2004')/100;
      Abs=xlsread('data_M', 'optical properties', 'E3:E2004')/100; 
    end
    
    if strcmp(glass,'nanofluid_1')==1
      Tr=xlsread('data_M', 'optical properties', 'R3:R2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'S3:S2004')/100;
      Abs=xlsread('data_M', 'optical properties', 'T3:T2004')/100; 
    end
    
      if strcmp(glass,'nanofluid_2')==1
      Tr=xlsread('data_M', 'optical properties', 'X3:X2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'Y3:Y2004')/100;
      Abs=xlsread('data_M', 'optical properties', 'Z3:Z2004')/100; 
      end
    
     if strcmp(glass,'nanofluid_3')==1
      Tr=xlsread('data_M', 'optical properties', 'AD3:AD2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'AE3:AE2004')/100;
      Abs=xlsread('data_M', 'optical properties', 'AF3:AF2004')/100; 
     end
  
     
     if strcmp(glass,'valvoline')==1
      Tr=xlsread('data_M', 'optical properties', 'AK3:AK2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'AL3:AL2004')/100;
      Abs=xlsread('data_M', 'optical properties', 'AM3:AM2004')/100; 
     end
     
     if strcmp(glass,'glycerol')==1
      Tr=xlsread('data_M', 'optical properties', 'AQ3:AQ2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'AR3:AR2004')/100;
      Abs=xlsread('data_M', 'optical properties', 'AS3:AS2004')/100; 
     end
     
         
     if strcmp(glass,'fake')==1
      Tr=xlsread('data_M', 'optical properties', 'AW3:AW2004')/100; 
      Rel=xlsread('data_M', 'optical properties', 'AX3:AX2004')/100;
      Abs=xlsread('data_M', 'optical properties', 'AY3:AY2004')/100; 
     end
     
    end
