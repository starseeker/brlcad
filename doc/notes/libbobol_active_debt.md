# libBObol active debt

Last reviewed: 2026-07-18

This is the complete active-debt list for the direct libged/libBObol drawing
contract.  Completed migration history belongs in the coverage notes and must
not be interpreted as an alternate supported path.

| Owner | Work | Enabling or removal condition | Gate |
|---|---|---|---|
| Platform maintainers | Windows core/headless and WGL/Tk/Qt validation | Native CI with package, endpoint replacement, resize/input/capture, and teardown coverage | Platform matrix moves the rows from deferred only after native runs pass |
| Platform maintainers | macOS core/headless and native host | Choose and support a Coin/Obol context provider, then add sanitizer and lifecycle CI | Platform matrix moves the rows from deferred only after native runs pass |
| libBObol host maintainers | Direct custom `BObolWindowHost` binding lacks a factory capability descriptor | Remove after supported external custom hosts register `bobol_host_factory` implementations | `libBObol_window_host` covers ownership and policy restrictions meanwhile |
| Rendering maintainers | Stereo presentation is not advertised by a validated host | Validate hardware capability, mode-specific output, capture, and mono restoration; software hosts must reject it until separately proven | Add host-specific image/lifecycle coverage before exposing policy |
| ART workstream | Appleseed tile publication is source-audited but not runtime-validated in this build | Run the imgstream publication path in an Appleseed-enabled build | Add the runtime result to the legacy-dependency inventory |
| LoD maintainers | FAA Generic_Twin full-detail refinement is too expensive for the routine suite | Keep refinement within the bounded progressive provider model | Promote the opt-in real-model workflow to a default gate |
| CI maintainers | TSan and LSan cannot start under some ptrace/restricted containers | Run them on a compatible native Linux worker; restricted workers validate the complete instrumented build only | Focused `bobol_sanitizer` runtime must pass on the compatible worker |

No branch-local controller attachment adapter, evaluated-wire compatibility
layer, alternate scene owner, or legacy redraw fallback is deferred here.
Reintroducing one requires an architecture decision, not an entry in this list.
