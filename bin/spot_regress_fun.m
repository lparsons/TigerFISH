function  [sl Mean Rsq] =  spot_regress_fun( p, dim)
clear x_7
persistent x_7 y_7 z_7    
persistent x_9 y_9 z_9
persistent x_11 y_11 z_11
persistent A7  A9 A11
if  isempty(x_7) 
    [x_7 y_7 z_7] = meshgrid( 1:7 );  A7 = ones( 7^3, 2 );
    [x_9 y_9 z_9] = meshgrid( 1:9 );  A9 = ones( 9^3, 2 );
    [x_11 y_11 z_11] = meshgrid( 1:11, 1:11, 1:7 );   A11 = ones( 7*11^2, 2 );  

    x_7=x_7(:); y_7=y_7(:); z_7=z_7(:);
    x_9=x_9(:); y_9=y_9(:); z_9=z_9(:);
    
    x_11=x_11(:); y_11=y_11(:); z_11=z_11(:);
end 


p = p(:) * ( 1/sum( p(:) ) );    

switch dim
    case 7       
    Mean(1) = sum( p .*  x_7 ); 
    Mean(2) = sum( p .*  y_7 );
    Mean(3) = sum( p .*  z_7 );

    x_gr7 = x_7 - Mean(1);
    y_gr7 = y_7 - Mean(2);
    z_gr7 = z_7 - Mean(3);

    A7(:,1) = sqrt( x_gr7.^2 + y_gr7.^2 + z_gr7.^2 );  %p=log(p);
    sl = linsolve( A7, p );
    hat_p = A7*sl;

    case 9     
    Mean(1) = sum( p .*  x_9 ); 
    Mean(2) = sum( p .*  y_9 );
    Mean(3) = sum( p .*  z_9 );

    x_gr9 = x_9 - Mean(1);
    y_gr9 = y_9 - Mean(2);
    z_gr9 = z_9 - Mean(3);            


    A9(:,1) = sqrt( x_gr9.^2 + y_gr9.^2 + z_gr9.^2 );  %p=log(p);
    sl = linsolve( A9, p );
    hat_p = A9*sl;
    
    case 11     
    Mean(1) = sum( p .*  x_11 ); 
    Mean(2) = sum( p .*  y_11 );
    Mean(3) = sum( p .*  z_11 );

    x_gr11 = x_11 - Mean(1);
    y_gr11 = y_11 - Mean(2);
    z_gr11 = z_11 - Mean(3);            


    A11(:,1) = sqrt( x_gr11.^2 + y_gr11.^2 + z_gr11.^2 );  %p=log(p);
    sl = linsolve( A11, p );
    hat_p = A11*sl;
end 

Rsq = 1-sum( (p-hat_p).^2 ) / ...
        sum( (p-mean(p)).^2 );     

%fprintf( 'Rsq: %1.3f \t %1.2e\t%1.2e\n', Rsq, sl );















