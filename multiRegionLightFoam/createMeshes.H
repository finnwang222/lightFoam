    regionProperties rp(runTime);
    
    PtrList<fvMesh> calRegions(rp.regionNames().size());

    forAll(rp.regionNames(), i)
    {
        Info<< "Create mesh for region " << rp.regionNames()[i]
            << " for time = " << runTime.timeName() << nl << endl;

        calRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    rp.regionNames()[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }


