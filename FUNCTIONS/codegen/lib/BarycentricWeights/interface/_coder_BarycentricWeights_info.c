/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_BarycentricWeights_info.c
 *
 * Code generation for function 'BarycentricWeights'
 *
 */

/* Include files */
#include "_coder_BarycentricWeights_info.h"

/* Function Definitions */
const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  const char * data[8] = {
    "789ced58cf6f1241149e36b5da8366a3d193d1c683079bb02d6d84f6547e2c8194d242694b699a74d91d60e2cc0ece2e089efa3f7835f1ea4d8f1efd13bce8c1"
    "937f8a4b97a5ec9809287469ebbe840c6fbfe57d6f5fde7c6f193093d99e0100dc038e7d7be2ac777bbed45b6781d7787ca6b7dee27cd0bf3ee7f99d8bbfefad",
    "1a352cd8b61c072303e69aa40299ed182a81fd303a25c8500dabd86940c0a049710beae748156158440466e9809346b643520350dfe942ddef893ad45eed3509"
    "6075f3225d3ce88081fa7c143cffdc88f54908ea2371f8b17292dc900bb0414d645186a079e1745294ed2615534eede712c5cc4e6e4f8eabaca341c362483b84",
    "a856b7cc10f1e47d2ac86b7ec4bcf9d5b5057067c07bbce9f2b505f146add343019fc4e11ad5210b21bb7998a1e2108646cdaa83e1cf3d2c0fde4479b8e6f27d"
    "fe473e377e6e089f8b1f67b225e524b121ef325a632a59ecf6b2296fc78ad9585c2e8497572215d9a21457685b86049f7f96ceab252fb9e5929d72755bc5c73e",
    "497d7af023e66f5ffabd0fa6c7e7d7be2b660fcaeb442fb1d8e15aa4bdb27a4096735b0379ec0ee119960710f87ec50f74de6bde7e0b6f8eabaf6efcdb023ea9"
    "87301db5900ec7d67337fe3ce75ff039884e9b150c27a7e719219f17ff6b3dc7a82213d5c26a45a60d53eed5a9d70d7ef6c3d9bbf9fcf740cf2f89cf2f3d8fbe",
    "c957f74b6b2f23e9162eac27a3d12cddd7c0cdd1f3ebbe8f1b82fc46edbb59ce776d81bbdfb1679bcefabcaff36782f8a3d6efa9805fe270ee7d1e99f126c256"
    "c6b0ff0e427be84d6d0e7c1993af24e4f3e29379afe7cb16223ece830fbf7e06f3e0aace8347023e89c3775aaf13b16aba152de6d73a9a428c70aa8cd3c13cb8",
    "2af3e05490df64fbeec5c4f47f51c0277138a7ffa6a662d86e242869a816b215795afaff754cbe23219f179f8cfeff51b669f44f70ce73897c7ecd81bade60e5"
    "d25b45ab4592d94867b5a9c6f2d1783007feaf39b031b173fdfb023e89c3b939603fbd73fdba9e036d09f9bcf864f4df2e97ffe78281de5f229f5f7a5f8bc73b",
    "7a645551caed702bbd5544c92356b801effdbf01f52fd00b", "" };

  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(data, 7912U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xResult;
  mxArray *xEntryPoints;
  const char * fldNames[4] = { "Name", "NumberOfInputs", "NumberOfOutputs",
    "ConstantInputs" };

  mxArray *xInputs;
  const char * b_fldNames[4] = { "Version", "ResolvedFunctions", "EntryPoints",
    "CoverageInfo" };

  xEntryPoints = emlrtCreateStructMatrix(1, 1, 4, fldNames);
  xInputs = emlrtCreateLogicalMatrix(1, 1);
  emlrtSetField(xEntryPoints, 0, "Name", emlrtMxCreateString(
    "BarycentricWeights"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs", emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs", emlrtMxCreateDoubleScalar
                (1.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  xResult = emlrtCreateStructMatrix(1, 1, 4, b_fldNames);
  emlrtSetField(xResult, 0, "Version", emlrtMxCreateString(
    "9.3.0.713579 (R2017b)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions", (mxArray *)
                emlrtMexFcnResolvedFunctionsInfo());
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_BarycentricWeights_info.c) */
