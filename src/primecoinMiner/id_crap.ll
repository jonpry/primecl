; ModuleID = 'get_global_size.ll'

declare i32 @llvm.r600.read.global.size.x() nounwind readnone

declare i32 @llvm.r600.read.global.size.y() nounwind readnone

declare i32 @llvm.r600.read.global.size.z() nounwind readnone

define i32 @get_global_size(i32 %dim) nounwind readnone alwaysinline {
  switch i32 %dim, label %default [
    i32 0, label %x_dim
    i32 1, label %y_dim
    i32 2, label %z_dim
  ]

x_dim:                                            ; preds = %0
  %x = call i32 @llvm.r600.read.global.size.x() nounwind readnone
  ret i32 %x

y_dim:                                            ; preds = %0
  %y = call i32 @llvm.r600.read.global.size.y() nounwind readnone
  ret i32 %y

z_dim:                                            ; preds = %0
  %z = call i32 @llvm.r600.read.global.size.z() nounwind readnone
  ret i32 %z

default:                                          ; preds = %0
  ret i32 0
}

declare i32 @llvm.r600.read.tgid.x() nounwind readnone

declare i32 @llvm.r600.read.tgid.y() nounwind readnone

declare i32 @llvm.r600.read.tgid.z() nounwind readnone

define i32 @get_group_id(i32 %dim) nounwind readnone alwaysinline {
  switch i32 %dim, label %default [
    i32 0, label %x_dim
    i32 1, label %y_dim
    i32 2, label %z_dim
  ]

x_dim:                                            ; preds = %0
  %x = call i32 @llvm.r600.read.tgid.x() nounwind readnone
  ret i32 %x

y_dim:                                            ; preds = %0
  %y = call i32 @llvm.r600.read.tgid.y() nounwind readnone
  ret i32 %y

z_dim:                                            ; preds = %0
  %z = call i32 @llvm.r600.read.tgid.z() nounwind readnone
  ret i32 %z

default:                                          ; preds = %0
  ret i32 0
}

declare i32 @llvm.r600.read.tidig.x() nounwind readnone

declare i32 @llvm.r600.read.tidig.y() nounwind readnone

declare i32 @llvm.r600.read.tidig.z() nounwind readnone

define i32 @get_local_id(i32 %dim) nounwind readnone alwaysinline {
  switch i32 %dim, label %default [
    i32 0, label %x_dim
    i32 1, label %y_dim
    i32 2, label %z_dim
  ]

x_dim:                                            ; preds = %0
  %x = call i32 @llvm.r600.read.tidig.x() nounwind readnone
  ret i32 %x

y_dim:                                            ; preds = %0
  %y = call i32 @llvm.r600.read.tidig.y() nounwind readnone
  ret i32 %y

z_dim:                                            ; preds = %0
  %z = call i32 @llvm.r600.read.tidig.z() nounwind readnone
  ret i32 %z

default:                                          ; preds = %0
  ret i32 0
}

declare i32 @llvm.r600.read.local.size.x() nounwind readnone

declare i32 @llvm.r600.read.local.size.y() nounwind readnone

declare i32 @llvm.r600.read.local.size.z() nounwind readnone

define i32 @get_local_size(i32 %dim) nounwind readnone alwaysinline {
  switch i32 %dim, label %default [
    i32 0, label %x_dim
    i32 1, label %y_dim
    i32 2, label %z_dim
  ]

x_dim:                                            ; preds = %0
  %x = call i32 @llvm.r600.read.local.size.x() nounwind readnone
  ret i32 %x

y_dim:                                            ; preds = %0
  %y = call i32 @llvm.r600.read.local.size.y() nounwind readnone
  ret i32 %y

z_dim:                                            ; preds = %0
  %z = call i32 @llvm.r600.read.local.size.z() nounwind readnone
  ret i32 %z

default:                                          ; preds = %0
  ret i32 0
}

declare i32 @llvm.r600.read.ngroups.x() nounwind readnone

declare i32 @llvm.r600.read.ngroups.y() nounwind readnone

declare i32 @llvm.r600.read.ngroups.z() nounwind readnone

define i32 @get_num_groups(i32 %dim) nounwind readnone alwaysinline {
  switch i32 %dim, label %default [
    i32 0, label %x_dim
    i32 1, label %y_dim
    i32 2, label %z_dim
  ]

x_dim:                                            ; preds = %0
  %x = call i32 @llvm.r600.read.ngroups.x() nounwind readnone
  ret i32 %x

y_dim:                                            ; preds = %0
  %y = call i32 @llvm.r600.read.ngroups.y() nounwind readnone
  ret i32 %y

z_dim:                                            ; preds = %0
  %z = call i32 @llvm.r600.read.ngroups.z() nounwind readnone
  ret i32 %z

default:                                          ; preds = %0
  ret i32 0
}
