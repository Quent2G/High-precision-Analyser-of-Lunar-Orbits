% metakernelcheck builds a metakernel.tm file to load kernels in the spice library.
% WARNING: such a file depends on ORBpropagation.m and it should not be run separately
% If you run only this file, a folder ORB will be created in the current directory
% where the file metakernel.tm will be saved. If you call cspice_furnsh('metakernel.tm')
% all the files listed in the "strend" variable, and that have to be included in a "ker"
% folder inside the "ORB" folder, will be loaded in the SPICE tool kernel.
% 
function [] = metakernelcheck()
    file = fopen('metakernel.tm','w');
    strfirst='KPL/MK \n \\begindata \n\t\t PATH_VALUES       = (\n\t\t\t''';
    filespath = regexprep([pwd,'/'],'\','/');
    strend=['ker''\n\t\t)\n\t\tPATH_SYMBOLS      = (''DATA'')\n\t\tKERNELS_TO_LOAD   = (',...
                    '\n\t\t\t''$DATA/de430.bsp'',', ...
                    '\n\t\t\t''$DATA/pck00010.tpc'',', ...
                    '\n\t\t\t''$DATA/naif0012.tls'',',...
                    '\n\t\t\t''$DATA/moon_pa_de421_1900-2050.bpc'',',...
                    '\n\t\t\t''$DATA/moon_080317.tf'',',...
                    '\n\t\t\t''$DATA/moon_assoc_me.tf'')',...
                    '\n\\begintext'];
    strdoc = [strfirst,filespath,strend];
    fprintf(file,strdoc);
    fclose(file);
end
    
