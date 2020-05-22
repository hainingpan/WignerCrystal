function [added] = mathOpsIntegrated(in1, in2,len)
%#codegen
% for code generation, preinitialize the output variable
% data type, size, and complexity 
added = 0;
% generate an include in the C code
coder.cinclude('adder.h');
% evaluate the C function
added = coder.ceval('adder', in1, in2,len); 
end