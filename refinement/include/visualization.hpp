#ifndef VISUALIZATION_HPP
#define VISUALIZATION_HPP

#include "refinementTree.hpp"

#include "output/graphics.hpp"

struct ZeroState
{
    template< typename Rtree, typename IterT >
    double operator ()( const Rtree& rtree, const IterT& cexBegin, const IterT& cexEnd, const IterT& istate ) const
    {
	return 0;
    }

    std::string name() const { return "zero_state";}
};

struct ZeroCex
{
    double operator ()( const double& a, const double& e ) const
    {
	return 0;
    }

    std::string name() const
    {
	return "zero_counterexample";
    }
};

template< typename Rtree, typename InisetT >
Ariadne::Colour stateColor( const Rtree& rtree, const InisetT& initialSet, const AbstractInteriorTreeValue< typename Rtree::EnclosureT >& tv, bool reachable )
{
    if( reachable )
    {
	if( definitely( tv.isSafe() ) )
	{
	    if( possibly( !initialSet.separated( tv.getEnclosure() ) ) )
		return Ariadne::cyan;
	    else
		return Ariadne::green;
	}
	else if( definitely( !tv.isSafe() ) )
	    return Ariadne::red;
	else
	    return Ariadne::Colour( 1, 0.7, 0 );
    }
    return Ariadne::white;
}

template< typename E, typename InisetT >
Ariadne::Figure visualize( const RefinementTree< E >& rtree, const InisetT& initialSet, const Ariadne::Effort& effort )
{
    auto ncomp = typename RefinementTree< E >::NodeComparator( rtree );
    VisitMap< typename RefinementTree< E >::EnclosureT > reachMap( ncomp );
    ZeroState stateH; ZeroCex cexH;
    CounterexampleStore< typename RefinementTree< E >::EnclosureT, ZeroState, ZeroCex > cstore( stateH, cexH );
    std::function< Ariadne::ValidatedUpperKleenean( const typename RefinementTree< E >::EnclosureT&, const Ariadne::BoundedConstraintSet& ) > interPred =
	[&effort] (auto& enc, auto& cset) {return !(cset.separated( enc ).check( effort ) ); };
    auto initialAbs = rtree.intersection( initialSet, interPred );
    findCounterexample( rtree, initialAbs.begin(), initialAbs.end(), cstore, reachMap );

    Ariadne::Figure fig( rtree.initialEnclosure(), Ariadne::PlanarProjectionMap( 2,0,1 ) );
    for( auto vs = graph::vertices( rtree.leafMapping() ); vs.first != vs.second; ++vs.first )
    {
	auto vval = rtree.nodeValue( *vs.first );
	if( vval )
	{
	    auto tv = vval.value().get();
	    auto imapEntry = reachMap.find( *vs.first );
	    fig.set_fill_colour( stateColor( rtree, initialSet, tv, imapEntry != reachMap.end() && imapEntry->second ) );
	    fig.draw( tv.getEnclosure() );
	}
    }
    return fig;
}

#endif
