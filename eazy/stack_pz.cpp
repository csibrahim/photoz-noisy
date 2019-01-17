#include <vif.hpp>

using namespace vif;

int vif_main(int argc, char* argv[]) {
    vec1s pz_file;
    vec1d z_spec;
    fits::read_table(argv[1], ftable(pz_file, z_spec));

    vec1d dzgrid = rgen_step(-9.0, 9.0, 0.001);
    vec1d pstack(dzgrid.size());

    auto pg = progress_start(pz_file.size());
    for (uint_t i : range(pz_file)) {
        vec1d z, pz;
        ascii::read_table(pz_file[i], z, _, _, pz);

        // Cleanup
        pz[where(pz < 0)] = 0;

        z -= z_spec[i];
        vec1d ipz = interpolate(pz, z, dzgrid);
        ipz[where(dzgrid < min(z) || dzgrid > max(z) || ipz < 0)] = 0;
        ipz /= total(ipz);

        pstack += ipz;

        progress(pg);
    }

    pstack /= total(pstack);

    fits::write_table(argv[2], ftable(dzgrid, pstack));

    return 0;
}
