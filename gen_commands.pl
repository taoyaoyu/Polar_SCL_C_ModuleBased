#!/usr/bin/perl

# Open new file to write
open(commands, ">commands.txt");

my @snr = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5);
my @num_err = (10000, 10000, 5000, 1000, 200, 100, 50, 20);
my @l_value = (1,2,4);
my @integer_bits = (5, 6, 8);
my @fractional_bits = (0, 0, 0);

printf commands "gcc -std=c99 -g -O3 -o polar_sim main.c src/*.c -lm\n\n";

for(my $i = 0; $i < 3; $i++){
	for(my $k = 0; $k < 3; $k++){
		for (my $j = 0; $j < 8; $j++){
			printf commands "polar_sim %d %1.1f %d %d %d\n", $num_err[$j], $snr[$j], $l_value[$i], $integer_bits[$k], $fractional_bits[$k];
		}
		printf commands "\n";
	}
	printf commands "\n";
}

close( commands );
