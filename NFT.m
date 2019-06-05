function varargout = NFT(varargin)

if nargout
    [varargout{1:nargout}] = Neuroelectromagnetic_Forward_Modeling_Toolbox(varargin{1:nargin});
else
    Neuroelectromagnetic_Forward_Modeling_Toolbox(varargin{1:nargin});
end
