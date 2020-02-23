clear ;
addpath('C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out')
Nspecs = 62800 ;
for spec = 1:Nspecs
    if mod(spec, 100) == 0
        disp(['spec ' num2str(spec) '!'])
    end
    PH_GER(spec) ; 
end

