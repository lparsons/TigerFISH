function [r p] = pearp ( m,n, N )

row=size( m,1 ); 

r = pear( m,n ); 

if nargin <3 || isempty(N)
    N=1e4;
end

if nargout == 2  
    Rpermuted = zeros( N, 1 );
    for i=1:N
      Rpermuted(i) = pear( m( randperm(row) ),...
                           n( randperm(row) ) );
    end
    if r > 0
        p = sum( Rpermuted > r )/N;
    else
        p = sum( Rpermuted < r )/N;
    end
end

 




function r = pear( m, n )


[row clm] = size( m );  
            
for i=1:clm,                     
        m(:,i) = m(:,i) -  mean( m(:,i) );                                               
        m(:,i) = m(:,i) /   std( m(:,i) );   
        
        n(:,i) = n(:,i) -  mean( n(:,i) );                                               
        n(:,i) = n(:,i) /   std( n(:,i) ); 
end


r   =   m' * n / (row-1);

