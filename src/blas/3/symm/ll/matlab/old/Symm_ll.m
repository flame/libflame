function [ X ] = Symm_ll( variant, alpha, A, B, C, nb )
%
%
if ( nb == 1 )
  if ( variant == 1 ) 
    X = FLA_Symm_ll_unb_var1( alpha, A, B, C );
  elseif ( variant == 2 )
    X = FLA_Symm_ll_unb_var2( alpha, A, B, C );
  elseif ( variant == 3 )
    X = FLA_Symm_ll_unb_var3( alpha, A, B, C );
  elseif ( variant == 4 )
    X = FLA_Symm_ll_unb_var4( alpha, A, B, C );
  elseif ( variant == 5 )
    X = FLA_Symm_ll_unb_var5( alpha, A, B, C );
  elseif ( variant == 6 )
    X = FLA_Symm_ll_unb_var6( alpha, A, B, C );
  elseif ( variant == 7 )
    X = FLA_Symm_ll_unb_var7( alpha, A, B, C );
  elseif ( variant == 8 )
    X = FLA_Symm_ll_unb_var8( alpha, A, B, C );
  elseif ( variant == 9 )
    X = FLA_Symm_ll_unb_var9( alpha, A, B, C );
  elseif ( variant == 10 )
    X = FLA_Symm_ll_unb_var10( alpha, A, B, C );
  end
  else
  if ( variant == 1 ) 
    X = FLA_Symm_ll_blk_var1( alpha, A, B, C, nb );
  elseif ( variant == 2 )
    X = FLA_Symm_ll_blk_var2( alpha, A, B, C, nb );
  elseif ( variant == 3 )
    X = FLA_Symm_ll_blk_var3( alpha, A, B, C, nb );
  elseif ( variant == 4 )
    X = FLA_Symm_ll_blk_var4( alpha, A, B, C, nb );
  elseif ( variant == 5 )
    X = FLA_Symm_ll_blk_var5( alpha, A, B, C, nb );
  elseif ( variant == 6 )
    X = FLA_Symm_ll_blk_var6( alpha, A, B, C, nb );
  elseif ( variant == 7 )
    X = FLA_Symm_ll_blk_var7( alpha, A, B, C, nb );
  elseif ( variant == 8 )
    X = FLA_Symm_ll_blk_var8( alpha, A, B, C, nb );
  elseif ( variant == 9 )
    X = FLA_Symm_ll_blk_var9( alpha, A, B, C, nb );
  elseif ( variant == 10 )
    X = FLA_Symm_ll_blk_var10( alpha, A, B, C, nb );
  end
end

return;
