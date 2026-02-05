from SummandTerm import MultipoleInteraction



s = MultipoleInteraction(multipole_order_1 = 1, multipole_order_2 = 1)
print( s.full_string_tensor() )
s.simplify()
print( " <==> \t", s.full_string_tensor() )

print()


s = MultipoleInteraction(multipole_order_1 = 1, multipole_order_2 = 2)
print( s.full_string_tensor() )
s.simplify()
print( " <==> \t", s.full_string_tensor() )