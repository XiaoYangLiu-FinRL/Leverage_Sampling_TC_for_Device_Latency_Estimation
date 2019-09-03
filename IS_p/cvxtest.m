    cvx_begin

    variable x(n1,n2)

    minimize( norm( L * x,'fro' ))

    subject to

        

    cvx_end