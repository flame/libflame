nvariants = 10
m = 20
n = 10
nb = 4

% create a random matrix A

A = rand( m, m );

% create a random matrix B 

B = rand( m, n );

% create a random matrix C 

C = rand( m, n );

% Create a table to report results

table_ll = zeros( nvariants, 2 );

% Compute D = symm( A ) * B + C with matlab

D = ( tril( A ) + tril( A, -1 )' ) * B + C;

for variant=1:nvariants
% Compute E = symm( A ) * B + C with unblocked variant

  E = Symm_ll( variant, 1.0, A, B, C, 1 );

  table_ll( variant, 1 ) = norm( D - E, 1 );

% Compute E = symm( A ) * B + C with blocked variant

  E = Symm_ll( variant, 1.0, A, B, C, nb );

  table_ll( variant, 2 ) = norm( D - E, 1 );
end

printf(" Results for Symm_ll \n\n");
printf("        |      Difference       \n");
printf("        +--------------------------\n");
printf("variant |  unblocked  |   blocked   \n");
printf("========+=============+=============\n");
for variant=1:nvariants
  printf("  %2d    | %5.2e    |  %5.2e\n", ...
                           variant, ...
                           table_ll( variant, 1 ), ...
                           table_ll( variant, 2 ) );
end


  
