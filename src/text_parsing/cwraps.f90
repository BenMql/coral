Module cwraps
      Use, intrinsic ::  Iso_C_Binding
      Interface
         Integer Function parse_text(compt, tot, l, aol) bind(C, name="parse_text")
               Use, intrinsic ::  Iso_C_Binding
               Integer :: compt, tot
               Type(C_ptr), Intent(In), Value :: l
               Type(C_ptr), Intent(In), Value :: aol
         End Function parse_text
         Integer Function parse_output(compt, tot, l, aol) bind(C, name="parse_output")
               Use, intrinsic ::  Iso_C_Binding
               Integer :: compt, tot
               Type(C_ptr), Intent(In), Value :: l
               Type(C_ptr), Intent(In), Value :: aol
         End Function parse_output
         Integer Function parse_timeseries(compt, tot, l, aol) bind(C, name="parse_timeseries")
               Use, intrinsic ::  Iso_C_Binding
               Integer :: compt, tot
               Type(C_ptr), Intent(In), Value :: l
               Type(C_ptr), Intent(In), Value :: aol
         End Function parse_timeseries
      End Interface
End Module cwraps
