function outputflag=FSP2SB3(system, name, theta, p0)
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
       Str=[Str 'd/dt(p' num2str(i) ') = '];
       for j=1:size(A,2)
           if (~isempty(A{i,j}))
           Str = [Str '+(' A{i,j} ')*p' num2str(j)];
           end
       end
           Str=[Str ' \n '];
   end
   toc
   for i=1:size(A,1)
       Str=[Str ' \n p' num2str(i) '(0)=' num2str(p0(i))];
   end
   output=['********** MODEL NAME \n ' name ' \n ********** MODEL NOTES \n '];
   output=[output '********** MODEL STATES \n ' Str ' \n '];
   output=[output '********** MODEL PARAMETERS \n '];
   for i=1:length(theta)
       output=[output 'theta' num2str(i) '=' num2str(theta(i)) ' \n '];
   end
   output=[output '********** MODEL VARIABLES \n '];
   output=[output '********** MODEL REACTIONS \n '];
   output=[output '********** MODEL FUNCTIONS \n '];
   output=[output '********** MODEL EVENTS \n '];
   output=[output '********** MODEL MATLAB FUNCTIONS'];
   f=fopen([name '.txt'], 'w');
   fprintf(f,output);
   fclose(f);
   outputflag=1;
   toc
end
