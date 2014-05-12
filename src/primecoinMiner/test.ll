; ModuleID = 'foo.cl'
target datalayout = "e-p:32:32-p1:64:64-p2:64:64-p3:32:32-p4:32:32-p5:64:64-i64:64-v16:16-v24:32-v32:32-v48:64-v96:128-v192:256-v256:256-v512:512-v1024:1024-v2048:2048-n32"
target triple = "r600--"

; Function Attrs: nounwind
define void @if_eq(i32 addrspace(1)* nocapture %out, i32 %arg0, i32 %arg1) #0 {
entry:
  %cmp = icmp eq i32 %arg0, %arg1
  %cond = zext i1 %cmp to i32
  store i32 %cond, i32 addrspace(1)* %out, align 4, !tbaa !2
  ret void
}

attributes #0 = { nounwind "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "unsafe-fp-math"="false" "use-soft-float"="false" }

!opencl.kernels = !{!0}
!llvm.ident = !{!1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1, !1}

!0 = metadata !{void (i32 addrspace(1)*, i32, i32)* @if_eq}
!1 = metadata !{metadata !"clang version 3.5 (https://github.com/llvm-mirror/clang.git a9ebe3817733d64547eee399d75b16421681b1af) (https://github.com/llvm-mirror/llvm.git 0a2572c1f1bd322d1517e15135033be88afc6cd7)"}
!2 = metadata !{metadata !3, metadata !3, i64 0}
!3 = metadata !{metadata !"int", metadata !4, i64 0}
!4 = metadata !{metadata !"omnipotent char", metadata !5, i64 0}
!5 = metadata !{metadata !"Simple C/C++ TBAA"}
