function basisobj =  create_product_basis(sbasisobj, tbasisobj)
%CREATE_PRODUCT_BASIS creates a basis object for bi-variate functions.
%  A product basis basis object is constructed from
%  two univariate bases for expanding a bivariate function.  
%  Each univariate basis is contained in a 'basis' object.  
%  The basis for the bivariate functional data object consists
%  of all possible pairs of the two sets of univariate basis functions.

%  This function is identical to CREATE_TP_BASIS.

%  Arguments
%  SBASISOBJ ... a functional data basis object for the first  argument s
%  TBASISOBJ ... a functional data basis object for the second argument t

%  Returns
%  BASISOBJ ... a functional data object

basisobj = create_TP_basis(sbasisobj, tbasisobj);

