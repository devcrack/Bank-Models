cc          = g++
cflags      = -Iheaders -c
output      = -o 01Hard_Spheere
in_files  = main.cpp file_manage.cpp hard_sphere_dynamics.cpp hard_sphere_mono.cpp static_variables.cpp system_variables.cpp
base_names  = $(basename $(in_files))
object_files = $(foreach str,$(base_names),$(str).o)


object_files: $(in_files)
	$(cc) $(in_files) $(cflags)

exe: object_files
	$(cc) $(object_files) $(output)
clean:
	rm ./*.o


