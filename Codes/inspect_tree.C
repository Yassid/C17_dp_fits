// Print the first TTree's name + branch list (and a couple of sample values)
// from an InterpSolver ROOT file, so the contaminant-fingerprint macro can be
// written against the real schema.
void inspect_tree(const char* path) {
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { printf("CANNOT OPEN %s\n", path); return; }
    TIter next(f->GetListOfKeys());
    TKey* k; TTree* t = nullptr;
    while ((k = (TKey*)next())) {
        if (TString(k->GetClassName()).Contains("Tree")) {
            t = (TTree*)f->Get(k->GetName());
            printf("TREE: %s   entries=%lld\n", k->GetName(), t->GetEntries());
            break;
        }
    }
    if (!t) { printf("NO TREE in %s\n", path); f->ls(); return; }
    printf("BRANCHES:\n");
    TObjArray* br = t->GetListOfBranches();
    for (int i = 0; i < br->GetEntries(); ++i) {
        TBranch* b = (TBranch*)br->At(i);
        TLeaf* lf = b->GetLeaf(b->GetName());
        printf("  %-28s %s\n", b->GetName(),
               lf ? lf->GetTypeName() : "?");
    }
}
