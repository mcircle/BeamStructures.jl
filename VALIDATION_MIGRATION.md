# Validation extraction (local, not merged)

Base: `master`, `3adf245ace11a8f985ae958ef8f37e64e734c6c7`.

The former `validation/` directory is now in the adjacent
`BeamStructuresValidations/scripts/legacy/` checkout. Its tests moved to
`BeamStructuresValidations/test/legacy/`. Core numerical and derivative tests
remain here. The topology smoke workflow was retained as a disabled migration
reference in the new package. No mechanical implementation was changed.

The old studies have historical targets and defaults. They do not implement the
new three-to-five-node study contract. Follow the new package README for status.
