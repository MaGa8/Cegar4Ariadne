
#include "geometry/box.hpp"
//#include "geometry/zonotope.hpp"
#include "geometry/polyhedron.hpp"
#include "geometry/box.decl.hpp"
#include "geometry/function_set.hpp"

// why is there no single module header to import? reason?
// contains constants, variables, expressions
#include "expression/expression.hpp"
#include "expression/space.hpp"
#include "expression/valuation.hpp"
// contains number types characterizing precision of number
#include "numeric/logical.hpp"
#include "numeric/module.hpp"
// used to produce functions
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/domain.hpp"

// present in rig numcs example
// #include "algebra/algebra.hpp"
// #include "function/formula.hpp"

#include <random>

int main()
{
    // evaluate a real function
    {
	Ariadne::RealVariable x( "x" );
	Ariadne::RealConstant c( "c", 2 );
	Ariadne::RealExpression e = x;
	Ariadne::RealValuation val( {x|c} );
	std::cout << e << " = " << Ariadne::evaluate( e, val ) << std::endl;
    }

    // evaluate effective function (on reals?) 
    {
	Ariadne::Space< Ariadne::Real > fspc;
	Ariadne::RealVariable x( "x" );
	fspc.insert( x );
	Ariadne::RealConstant c( "c", 2 );
	// does this even work for univariate functions?
	// can only obtain scalar multivariate ones from make_function (it seems)
	Ariadne::EffectiveScalarFunction f = Ariadne::make_function( fspc, x * x );
	
	// how can i evaluate effective univar functions? using real?
	Ariadne::RealVector veval = { 2 };
	// Ariadne::Box< 
	std::cout << f << " = " << f.evaluate( veval ) << std::endl;
    }

    // evaluate vector valued function
    {
	auto idty = Ariadne::EffectiveVectorFunction::identity( Ariadne::EuclideanDomain( 2 ) );
	// where is the subscript operator defined for functions?
	// what does it mean?
	// auto a = idty[ 0 ]; auto b = idty[ 1 ];
	Ariadne::RealVariable a( "a" ), b( "b" );
	Ariadne::Space< Ariadne::Real > fspc;
	fspc.insert( a ); fspc.insert( b );

	// memory access error here
	Ariadne::EffectiveVectorFunction g = Ariadne::make_function( fspc, { a + b, a - b } );

	Ariadne::RealVector v = {1,2};
	std::cout << g << " = " << g.evaluate( v ) << std::endl;
    }

    // evaluate some function at center of some box
    {
    	Ariadne::EffectiveBoxType bx1 = { {-2,-2}, {+2,+2} }
	    , bx2 = { {-1,-1}, {1,1} }
	    , bx3 = { {0,0}, {1,1} };
    	// what is the difference between center and midpoint?

	Ariadne::RealVariable x( "x" ); Ariadne::RealVariable y( "y" );
	Ariadne::Space< Ariadne::Real > fspace( { x, y } );
	
	// what does the type of x, y need to be to construct function?
	// Ariadne::Vector< Ariadne::FloatDPValue > c = bx.centre();
	Ariadne::Vector< Ariadne::Real > c = bx1.centre();
    	Ariadne::EffectiveVectorFunction f = Ariadne::make_function( fspace, { x + y, y + x } );

	// Ariadne::Vector< Ariadne::EffectiveIntervalType > cres = f.evaluate( c );
	// real is effective type in Ariadne
	Ariadne::Vector< Ariadne::Real > cres = f.evaluate( c );
	std::cout << f << " = " << c << std::endl;
	std::cout << "center of box 1 in box 2? "
		  << bx2.contains( Ariadne::Point< Ariadne::Real >( c ) )
		  << std::endl;
	std::cout << "mapped center of box 1 in box 2? "
		  << bx2.contains( Ariadne::Point< Ariadne::Real >( cres ) )
		  << std::endl;
    	// center is on boundary -- what would happen if???
	std::cout << "mapped center of box 1 in box 3? "
		  << bx3.contains( Ariadne::Point< Ariadne::Real >( cres ) )
		  << std::endl;

	Ariadne::EffectiveIntervalType bx1I0 = bx1[ 0 ];
	// ??? Ariadne::UpperIntervalType bx1I0Upper( bx1I0 );
	// ??? Ariadne::Bounds< Ariadne::FloatDP > bx1I0Bounds( bx1I0.lower() );
	Ariadne::ExactBoxType exBox = { {-1,1}, {-1,2} };
	Ariadne::ExactIntervalType exBoxI0 = exBox[ 0 ];
	Ariadne::UpperInterval< Ariadne::FloatDP > exBoxI0Upper( exBoxI0 );
	Ariadne::UpperBoxType ubox( exBox.array() );
	Ariadne::UpperBoxType imagedBox = Ariadne::image( ubox, f );
	
    	// // no need to obtain box from result of function
    	// // because of effective paradigm
    	// // will be able to tell whether is inside any box if not on boundary
    }

    {
	Ariadne::RealVariable x( "x" ), y( "y" );
	Ariadne::Space< Ariadne::Real > vspace = {x, y};
	Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x,y} );

	auto a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 );
	auto b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
	auto expr = a + b;
	Ariadne::EffectiveConstraint c = ( expr <= 1 );
	Ariadne::EffectiveConstraintSet cs = {c};

	Ariadne::ExactBoxType bx1 = { {0,0}, {2,2} }, bx2 = { {1,1}, {2,2} };
	std::cout << "does " << bx1 << " intersect " << bx2 << "? "
		  << Ariadne::possibly( Ariadne::intersection( bx1, bx2 ).is_empty() ) << std::endl;

	Ariadne::Effort effort( 10 );
	Ariadne::ValidatedLowerKleenean doesCover = cs.covers( bx1 ).check( effort );
	// Ariadne::definitely( doesCover );
	definitely( cs.covers( bx1 ) );
	// std::cout << "does " << cs << " cover " << bx1 << "? "
	// 	  << Ariadne::definitely< Ariadne::ValidatedLowerKleenean >( cs.covers( bx1 ).check( effort ) ) << std::endl;
	// why does this not work? not convertible to bool? sure...
    }

    {
	std::uniform_real_distribution<> valDist( -1000, 1000 );
	std::default_random_engine mRandom = std::default_random_engine( std::random_device()() );

	Ariadne::RealConstant w( "w", Ariadne::Real( valDist( mRandom ) ) )
	    , h( "h", Ariadne::Real( valDist( mRandom ) ) );
	Ariadne::EffectiveScalarFunction x = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	    , y = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
	

	Ariadne::ConstraintSet cs( { x >= -w, x <= w, y >= -h, y <= h } );

	for( uint i = 0; i < 1000; ++i )
	{
	    Ariadne::ExactPoint pt( { valDist( mRandom ), valDist( mRandom ) } );
	    Ariadne::ExactBoxType bx( pt );
	    if( pt[ 0 ].raw() >= -w.get_d() && pt[ 0 ].raw() <= w.get_d() && pt[ 1 ].raw() >= -h.get_d() && pt[ 1 ].raw() <= h.get_d() &&
		possibly( !cs.covers( Ariadne::ExactBoxType( pt ) ) ) )
		std::cout << "error" << std::endl;
	}

	std::cout << "everything went fine" << std::endl;
	std::cout << "constraint bb " << cs.constraint_bounds() << std::endl << std::endl;
    }

    {
	std::uniform_real_distribution<> valDist( 0, 1000 );
	std::default_random_engine mRandom = std::default_random_engine( std::random_device()() );

	double wid = valDist( mRandom ), hig = valDist( mRandom );

	Ariadne::RealConstant w( "w", Ariadne::Real( wid ) )
	    , h( "h", Ariadne::Real( hig ) );
	Ariadne::EffectiveScalarFunction x = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	    , y = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
	
	Ariadne::ExactBoxType smallBox = { {-wid/2, wid/2}, {-hig/2, hig/2} }
	    , bigBox = { {-wid*2, wid*2}, {-hig*2, hig*2} };
	
	Ariadne::ConstraintSet cs( { -w <= x <= w, -h <= y <= h } );
	Ariadne::Effort e( 10 );

	for( auto constraint : cs.constraints() )
	{
	    Ariadne::ConstraintSet singleConstraint( constraint );
	    if( !definitely( singleConstraint.covers( smallBox ).check( e ) ) )
		std::cout << singleConstraint << " does not cover small box " << smallBox << std::endl;
	    if( !definitely( singleConstraint.covers( bigBox ).check( e ) ) )
		std::cout << singleConstraint << " does not cover big box " << bigBox << std::endl;
	}
    }

    // what is the logic type of comparing upper box type with exact box?
    {
	Ariadne::ExactBoxType ebx = { {-1, 1}, {-1, 1} };
	Ariadne::RealBox rbx( ebx );
	Ariadne::BoundedConstraintSet bset( rbx  );
	Ariadne::UpperBoxType ubx = bset.bounding_box();
	Ariadne::ValidatedKleenean equal = ubx == ubx;
    }

    // do basic ops with zonotopes
    {
	Ariadne::ExactBoxType zbx = { {-1, 1}, {-1, 1} };
	// Ariadne::Zonotope z( zbx );
	Ariadne::Polyhedron h( zbx );

	Ariadne::RealVariable x( "x" ), y( "y" );
	Ariadne::Space< Ariadne::Real > vspace = {x, y};
	Ariadne::ValidatedVectorFunction f = Ariadne::make_function( vspace, {x,y} );
	// Ariadne::image( z, f );
	
    }

    
    return 0;
}
						    
