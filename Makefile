CC = gcc
CFLAGS = -O3 -march=native -Wall -Wextra

shock_tube: shock_tube.c tables.h
	$(CC) $(CFLAGS) -o $@ shock_tube.c -lm

clean:
	rm -f shock_tube

tables.h: export_tables.m methane_U.m methane_Cv.m methane_soundspeed.m write_c_table.m
	octave --no-gui --eval "run('export_tables.m')"

.PHONY: clean
