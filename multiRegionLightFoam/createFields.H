    // Initialise field pointer lists
    
    PtrList<volScalarField> I_cal(calRegions.size());
    
    PtrList<light::lightModel> lightModels(calRegions.size());

    // Populate solid field pointer lists
    forAll(calRegions, i)
    {
        Info<< "*** Reading  properties for region "
            << calRegions[i].name() << nl << endl;

        I_cal.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "IDefault",
                    runTime.timeName(),
                    calRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                calRegions[i]
            )
        );

        lightModels.set
        (
            i,
            light::lightModel::New(I_cal[i])
        );
    }
    
    
    
