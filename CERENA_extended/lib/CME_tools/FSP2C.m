function outputflag=FSP2C(system, name, theta, p0)
tic
for i=1:length(system.sigma_A_i)
    for j=1:length(system.sigma_A_i)
        if(~isempty(strfind(char(system.sigma_A_i{i}), ['theta(' num2str(j) ')'])))
    s(i,j)=j;
        else 
            s(i,j)=0;
        end
    end
end
    
    for i=1:length(s)
        str{i}=['theta' num2str(s(i,(s(i,:)~=0)))];
    end
    A=cell(size(system.A_i{1},1), size(system.A_i{1},2));
    
        
        for i=1:size(system.A_i{1},1)
            for j=1:size(system.A_i{1},2)
                for k=1:length(str)
                    if(system.A_i{k}(i,j)~=0)
            A{i,j}= ...
                [A{i,j} ...
                '+' str{k} '*' '(' ... 
                num2str(system.A_i{k}(i,j)) ')'];
                    end
                end
            end
         end
   Str=[];
   for i=1:size(A,1)
       Str=[Str 'DDTvector[' num2str(i-1) '] = '];
       for j=1:size(A,2)
           if (~isempty(A{i,j}))
           Str = [Str '+(' A{i,j} ')*stateVector[' num2str(j-1) ']'];
           end
       end
           Str=[Str '; \n '];
   end
   toc
   StrIC=['double defaultICs[' num2str(size(A,1)) '] = { \n '];
   for i=1:(size(A,1)-1)
       StrIC=[StrIC num2str(p0(i)) ', \n '];
   end
   StrIC=[StrIC num2str(p0(i+1)) '};'];
   StrStateNames=['char *stateNames[2121] = { '];
   for i=1:(size(A,1)-1)
   StrStateNames=[StrStateNames ' \n "p' num2str(i) '",'];
   end
   StrStateNames=[StrStateNames '"p' num2str(i+1) '"};'];
   StrPar=['double defaultParam[' num2str(length(theta)) '] = {'];
   for i=1:(length(theta)-1)
       StrPar=[StrPar num2str(theta(i)) ','];
   end
   StrPar=[StrPar num2str(theta(i+1)) '};'];
   StrParN=['char *parameterNames[' num2str(length(theta)) '] = {'];
   for i=1:(length(theta)-1)
       StrParN=[StrParN '\n "theta' num2str(i) '",'];
   end
   StrParN=[StrParN '\n "theta' num2str(i+1) '"};'];
   outputh=['#include "mex.h" \n const int NRSTATES = ' num2str(size(A,1)) ';'];
   outputh=[outputh '\n const int NRPARAMETERS = ' num2str(length(theta)) ';']; 
   outputh=[outputh '\n const int NRVARIABLES = 0; \n const int NRREACTIONS = 0; \n const int NREVENTS = 0; \n'];
   outputh=[outputh '\n' StrIC ' \n ' StrPar ' \n ' StrStateNames ' \n ' StrParN ' \n '];
   outputh=[outputh '\n char *variableNames[1]; \n char *reactionNames[1]; \n char *eventNames[1];'];
   outputh=[outputh ' \n void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);'];
   outputh=[outputh ' \n void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);'];
   outputh=[outputh ' \n void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) \n { \n CVODEmex25(nlhs, plhs, nrhs, prhs); \n } \n'];
   
   StrParNa='double ';
   for i=1:(length(theta)-1)
   StrParNa=[StrParNa 'theta' num2str(i) ','];
   end
   StrParNa=[StrParNa 'theta' num2str(i+1) ';'];
   StrP='\n ';
   for i=1:(length(theta))
   StrP=[StrP 'theta' num2str(i) ' = paramdataPtr->parametervector[' num2str(i-1) ']; \n'];
   end
   output=[];
   output=[output '#include "stddef.h" \n #include "stdarg.h" \n '];
   output=[output '#include "math.h" \n #include "CVODEmex25.h" \n'];
   output=[output '#include "' name '.h" \n #include "splineSB.h" \n'];
   output=[output '#include "mexmathaddon.h" \n #include "kineticformulas.h" \n'];
   output=[output 'double time; \n'];
   output=[output 'void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector) \n '];
   output=[output '{ \n ' StrParNa  ' \n time = time_local; \n '];
   output=[output StrP ' \n if (DOflag == DOFLAG_DDT) { \n ' Str];
   output=[output ' \n } else if (DOflag == DOFLAG_VARREAC) { \n '];
   output=[output '} else if (DOflag == DOFLAG_EVENTS) { \n } else if (DOflag == DOFLAG_EVENTASSIGN) { \n } \n}'];
   f=fopen([name '.h'], 'w');
   fc=fopen([name '.c'], 'w');
   fprintf(f,outputh);
   fclose(f);
   fprintf(fc, output);
   fclose(fc);
   outputflag=1;
   toc
end
