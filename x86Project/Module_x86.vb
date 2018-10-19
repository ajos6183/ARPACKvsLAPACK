Module Module_x86

    Dim Timer As New Stopwatch

    Sub Main()

        Console.WriteLine("Solver".PadRight(40) + "Time (Seconds)")
        Console.WriteLine("------".PadRight(40) + "--------------")

        Dim Kg = CSparse.IO.MatrixMarketReader.ReadMatrix(Of Double)(System.IO.Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Kg"))
        Dim Ke = CSparse.IO.MatrixMarketReader.ReadMatrix(Of Double)(System.IO.Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Ke"))
        Ke.Keep(Function(i, j, aij) i <= j)
        Kg.Keep(Function(i, j, aij) i <= j)

        Timer.Reset()
        Timer.Start()
        For k = 0 To 100
            Dim n As Integer = Ke.RowCount
            Dim Ke_UpperTriangle = New Double(n * n - 1) {}
            Dim Kg_UpperTriangle_EigenVectors = New Double(n * n - 1) {}
            For r As Integer = 0 To n - 1
                For c As Integer = r To n - 1
                    Ke_UpperTriangle(r + c * n) = Ke.At(r, c)
                    Kg_UpperTriangle_EigenVectors(r + c * n) = Kg.At(r, c)
                Next
            Next
            Dim ggg As Double() = EigSymmetric1(Kg_UpperTriangle_EigenVectors, Ke_UpperTriangle)
            'If k = 0 Then Console.WriteLine(String.Join(", ", ggg.Select(Function(x) x.ToString("F5")).ToArray()))
        Next
        Timer.Stop()
        Console.WriteLine("LAPACK/dsygv_:".PadRight(40) + (Timer.ElapsedMilliseconds / 1000).ToString("F3"))

        Timer.Reset()
        Timer.Start()
        Dim TaskArray1 = New Task(100) {}
        For k = 0 To 100
            Dim Index = k
            TaskArray1(Index) = Task.Factory.StartNew(Sub()
                                                          Dim n As Integer = Ke.RowCount
                                                          Dim Ke_UpperTriangle = New Double(n * n - 1) {}
                                                          Dim Kg_UpperTriangle_EigenVectors = New Double(n * n - 1) {}
                                                          For r As Integer = 0 To n - 1
                                                              For c As Integer = r To n - 1
                                                                  Ke_UpperTriangle(r + c * n) = Ke.At(r, c)
                                                                  Kg_UpperTriangle_EigenVectors(r + c * n) = Kg.At(r, c)
                                                              Next
                                                          Next
                                                          Dim ggg As Double() = EigSymmetric1(Kg_UpperTriangle_EigenVectors, Ke_UpperTriangle)
                                                          'If k = 0 Then Console.WriteLine(String.Join(", ", ggg.Select(Function(x) x.ToString("F5")).ToArray()))
                                                      End Sub)
        Next
        Task.WaitAll(TaskArray1)
        Timer.Stop()
        Console.WriteLine("LAPACK/dsygv_ [Multi-threaded]:".PadRight(40) + (Timer.ElapsedMilliseconds / 1000).ToString("F3"))

        Timer.Reset()
        Timer.Start()
        For k = 0 To 100
            Dim n As Integer = Ke.RowCount
            Dim Ke_UpperTriangle = New Double(n * n - 1) {}
            Dim Kg_UpperTriangle_EigenVectors = New Double(n * n - 1) {}
            For r As Integer = 0 To n - 1
                For c As Integer = r To n - 1
                    Ke_UpperTriangle(r + c * n) = Ke.At(r, c)
                    Kg_UpperTriangle_EigenVectors(r + c * n) = Kg.At(r, c)
                Next
            Next
            Dim ggg As Double() = EigSymmetric2(Kg_UpperTriangle_EigenVectors, Ke_UpperTriangle)
            'If k = 0 Then Console.WriteLine(String.Join(", ", ggg.Select(Function(x) x.ToString("F5")).ToArray()))
        Next
        Timer.Stop()
        Console.WriteLine("LAPACK/dsygvd_:".PadRight(40) + (Timer.ElapsedMilliseconds / 1000).ToString("F3"))

        Timer.Reset()
        Timer.Start()
        Dim TaskArray2 = New Task(100) {}
        For k = 0 To 100
            Dim Index = k
            TaskArray2(Index) = Task.Factory.StartNew(Sub()
                                                          Dim n As Integer = Ke.RowCount
                                                          Dim Ke_UpperTriangle = New Double(n * n - 1) {}
                                                          Dim Kg_UpperTriangle_EigenVectors = New Double(n * n - 1) {}
                                                          For r As Integer = 0 To n - 1
                                                              For c As Integer = r To n - 1
                                                                  Ke_UpperTriangle(r + c * n) = Ke.At(r, c)
                                                                  Kg_UpperTriangle_EigenVectors(r + c * n) = Kg.At(r, c)
                                                              Next
                                                          Next
                                                          Dim ggg As Double() = EigSymmetric2(Kg_UpperTriangle_EigenVectors, Ke_UpperTriangle)
                                                          'If k = 0 Then Console.WriteLine(String.Join(", ", ggg.Select(Function(x) x.ToString("F5")).ToArray()))
                                                      End Sub)
        Next
        Task.WaitAll(TaskArray2)
        Timer.Stop()
        Console.WriteLine("LAPACK/dsygvd_ [Multi-threaded]:".PadRight(40) + (Timer.ElapsedMilliseconds / 1000).ToString("F3"))

        Timer.Reset()
        Timer.Start()
        For k = 0 To 100
            Dim Eig As CSparse.Double.Solver.Arpack = New CSparse.Double.Solver.Arpack(Kg, Ke, True) With {.Tolerance = 10 ^ -5, .ComputeEigenVectors = True, .Iterations = 1000}
            Dim Solution As CSparse.Interop.ARPACK.ArpackResult(Of Double) = Eig.SolveGeneralized(Math.Min(Ke.ColumnCount - 2, 12), CSparse.Interop.ARPACK.Job.LargestMagnitude)
            'If k = 0 Then Console.WriteLine(String.Join(", ", Solution.EigenValues.Select(Function(x) x.Real.ToString("F5")).ToArray()))
        Next
        Timer.Stop()
        Console.WriteLine("ARPACK:".PadRight(40) + (Timer.ElapsedMilliseconds / 1000).ToString("F3"))
        Console.ReadKey()

    End Sub

    <Runtime.InteropServices.DllImport("liblapack", EntryPoint:="dsygv_", CallingConvention:=Runtime.InteropServices.CallingConvention.Cdecl), Security.SuppressUnmanagedCodeSecurity>
    Private Sub Lapack_dsygv(ByRef itype As Integer, ByRef jobz As Char, ByRef uplo As Char, ByRef n As Integer, ByVal A As Double(), ByRef lda As Integer,
                             ByVal B As Double(), ByRef ldb As Integer, ByVal w As Double(), ByVal work As Double(), ByRef lwork As Integer, ByRef info As Integer)
    End Sub
    Public Sub Dsygv(itype As Integer, jobz As Char, uplo As Char, n As Integer, A As Double(), lda As Integer, B As Double(), ldb As Integer, w As Double(), ByRef info As Integer)
        Dim lwork As Integer = -1
        Dim work As Double() = New Double(0) {0.0}
        Lapack_dsygv(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info)
        If info <> 0 OrElse work(0) <= 0.0 Then Return
        lwork = CInt(work(0))
        work = New Double(lwork - 1) {}
        Lapack_dsygv(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info)
    End Sub
    Public Function EigSymmetric1(A_UpperTriangle_EigenVectors As Double(), B_UpperTriangle As Double()) As Double()
        Dim n As Integer = CInt(Math.Sqrt(A_UpperTriangle_EigenVectors.GetLength(0)))
        Dim info As Integer = -1
        Dim w As Double() = New Double(n - 1) {}
        Dsygv(1, "N"c, "U"c, n, A_UpperTriangle_EigenVectors, n, B_UpperTriangle, n, w, info)
        If info = 0 Then
            Return w
        ElseIf info < 0 Then
            Throw New Exception("eigSymm: invalid parameter #" & (-info))
        Else
            If info <= n Then
                Throw New Exception(String.Format("eigSymm: did not converge! {0} off-diagonal elements unequal 0", info))
            ElseIf info < 2 * n Then
                Throw New Exception("eigSymm: B must be positive definite!")
            Else
                Throw New Exception("eigSymm: unknown error")
            End If
        End If
    End Function

    <Runtime.InteropServices.DllImport("liblapack", EntryPoint:="dsygvd_", CallingConvention:=Runtime.InteropServices.CallingConvention.Cdecl), Security.SuppressUnmanagedCodeSecurity>
    Private Sub Lapack_dsygvd(ByRef itype As Integer, ByRef jobz As Char, ByRef uplo As Char, ByRef n As Integer, ByVal A As Double(), ByRef lda As Integer,
                         ByVal B As Double(), ByRef ldb As Integer, ByVal w As Double(), ByVal work As Double(), ByRef lwork As Integer,
                              ByVal iwork As Integer(), ByRef llwork As Integer, ByRef info As Integer)
    End Sub
    Public Sub Dsygvd(itype As Integer, jobz As Char, uplo As Char, n As Integer, A As Double(), lda As Integer, B As Double(), ldb As Integer, w As Double(), ByRef info As Integer)
        Dim lwork As Integer = -1
        Dim work As Double() = New Double(0) {0.0}
        Dim ilwork As Integer = If(jobz = "V"c, 3 + 5 * n, 1)
        Dim iwork As Integer() = New Integer(ilwork - 1) {}
        Lapack_dsygvd(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, iwork, ilwork, info)
        If info <> 0 OrElse work(0) <= 0.0 Then Return
        lwork = CInt(work(0))
        work = New Double(lwork - 1) {}
        Lapack_dsygvd(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, iwork, ilwork, info)
    End Sub
    Public Function EigSymmetric2(A_UpperTriangle_EigenVectors As Double(), B_UpperTriangle As Double()) As Double()
        Dim n As Integer = CInt(Math.Sqrt(A_UpperTriangle_EigenVectors.GetLength(0)))
        Dim info As Integer = -1
        Dim w As Double() = New Double(n - 1) {}
        Dsygvd(1, "N"c, "U"c, n, A_UpperTriangle_EigenVectors, n, B_UpperTriangle, n, w, info)
        If info = 0 Then
            Return w
        ElseIf info < 0 Then
            Throw New Exception("eigSymm: invalid parameter #" & (-info))
        Else
            If info <= n Then
                Throw New Exception(String.Format("eigSymm: did not converge! {0} off-diagonal elements unequal 0", info))
            ElseIf info < 2 * n Then
                Throw New Exception("eigSymm: B must be positive definite!")
            Else
                Throw New Exception("eigSymm: unknown error")
            End If
        End If
    End Function

End Module