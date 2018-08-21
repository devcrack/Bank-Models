cc          = g++
cflags      = -Iheaders -c
output      = -o blds
in_files    = file_manage.cpp menu_things.cpp main.cpp
bases_names = $(basename $(in_files))
object_files = $(foreach str,$(bases_names),$(str).o)

object_files: $(in_files)
	$(cc) $(in_files) $(cflags)

exe: object_files
	$(cc) $(object_files) $(output)
.clean:
	rm ./*.o



