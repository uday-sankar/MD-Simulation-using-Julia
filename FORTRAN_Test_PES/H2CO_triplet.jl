# Get_PES_test.jl - ALL POINTERS VERSION
const LIBCH2O = "/Users/upm5017/Documents/Work/FORTRAN_Julia/Working files/libch2o_pes.so" # location of the compiled (as a shared library) fortran code
const PES_PATH = "/Users/upm5017/Documents/Work/FORTRAN_Julia/Working files" # location of the current folder

function ch2o_pes(x::Matrix{Float64})
    @assert size(x) == (4, 3) "Input must be 4×3 (H1,H2,C,O × x,y,z)"
    
    # 1. ALL SCALARS → Ref{} (gfortran pass-by-reference)
    igrad_ref = Ref{Cint}(0)                    # igrad=0 (energy only)
    p_ref = Ref{Cdouble}(0.0)                   # Output energy
    
    g_flat = zeros(Cdouble, 12)                 # 1×4×3 flattened
    d_flat = zeros(Cdouble, 12)                 # 1×1×4×3 flattened
    
    ## Path
    path_str = rpad(PES_PATH, 1024, ' ')

    # 3. ALL POINTERS ccall (gfortran convention)
    ccall((:pes_, LIBCH2O), Cvoid,
          (Ptr{Cdouble},     # x(4,3) → flattened
           Ptr{Cint},        # igrad → Ref
           Cstring,          # path (string is already pointer)
           Ptr{Cdouble},     # p(1) → Ref
           Ptr{Cdouble},     # g(1,4,3) → flattened
           Ptr{Cdouble}),    # d(1,1,4,3) → flattened
           x,                # Array pointer ✓
          igrad_ref,         # Scalar pointer ✓
          path_str,          # String pointer ✓
          p_ref,             # Output pointer ✓
          g_flat,            # Array pointer ✓
          d_flat)            # Array pointer ✓
    return p_ref[]
end

# Test (exact same coords as standalone Fortran)
function test_ch2o()
    print("===== !!!! Testing FORTRAN ccall !!!! =====")
    coords = [
        -0.628  0.000  0.936;   # H1
        -0.628  0.000 -0.936;   # H2  
         0.000  0.000  0.000;   # C
         1.210  0.000  0.000    # O
    ]
    energy = ch2o_pes(coords)
    println("\n Func Call Suecessfull \n CH2O Energy: $energy Hartree")
    print(" If Energy = -2.4502807685963415 Hartree Function excecution sucessfull\n ===== !!!! Done !!!! =====")
end

# Run test
test_ch2o()
